package lex1dx

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths}
import scala.jdk.CollectionConverters.*

class AgreeTest extends munit.FunSuite:

  private def builtin: Array[Double] =
    Array.tabulate(800)(i => math.sin(0.31 * i) + 0.05 * (i % 17) + ((i % 5) - 2) * 0.2)

  private def loadF64(path: Path): Option[Array[Double]] =
    if !Files.isRegularFile(path) then None
    else
      val bytes = Files.readAllBytes(path)
      if bytes.length == 0 || bytes.length % 8 != 0 then None
      else
        val n = bytes.length / 8
        val bb = java.nio.ByteBuffer.wrap(bytes).order(java.nio.ByteOrder.LITTLE_ENDIAN)
        val a = new Array[Double](n)
        var i = 0
        while i < n do { a(i) = bb.getDouble(); i += 1 }
        Some(a)

  private def formatResult(name: String, mode: String, kreq: Int, r: Result): String =
    val bounds = r.bounds.mkString(",")
    val centers = r.means.map(m => f"$m%.17g").mkString(",")
    s"case=$name mode=$mode k=$kreq k_eff=${r.k} D=${r.D} Dprime=${r.Dprime} sse=${f"${r.sse}%.17g"}\nbounds=$bounds\ncenters=$centers\n"

  test("write Scala dumps for the fixed set; compare to C++ if agree_cpp.txt is present"):
    val dir = Option(System.getenv("LEX1DX_AGREE_DIR")).map(Paths.get(_)).getOrElse(Paths.get("."))
    Files.createDirectories(dir)
    val out = new StringBuilder
    def run(name: String, xs: Array[Double], k: Int, mode: Mode): Result =
      val r = Lex.cluster(xs, k, Options(mode = mode))
      val ms = if mode == Mode.Sig5 then "sig5" else "tuned"
      out ++= formatResult(name, ms, k, r)
      out ++= "---\n"
      r

    val xs = builtin
    for k <- List(4, 8, 16, 64); mode <- List(Mode.Tuned, Mode.Sig5) do
      val r = run("builtin", xs, k, mode)
      assertEquals(r.k, k)

    val data = Option(System.getenv("LEX1DX_DATA")).map(Paths.get(_))
    data.foreach { root =>
      for fam <- List("uniform", "normal", "lognormal") do
        loadF64(root.resolve(s"synth/${fam}_n65536.f64")).filter(_.length == 65536).foreach { arr =>
          for k <- List(64, 256, 512); mode <- List(Mode.Tuned, Mode.Sig5) do
            run(s"${fam}_n65536", arr, k, mode)
        }
      loadF64(root.resolve("audio_s0.raw.f64")).filter(_.length >= 65536).foreach { arr =>
        val slice = arr.take(65536)
        for k <- List(64, 256, 512); mode <- List(Mode.Tuned, Mode.Sig5) do
          run("audio_slice65536", slice, k, mode)
      }
    }

    Files.writeString(dir.resolve("agree_scala.txt"), out.toString, StandardCharsets.US_ASCII)

    val cppPath = dir.resolve("agree_cpp.txt")
    if Files.isRegularFile(cppPath) then
      val cpp = Files.readString(cppPath, StandardCharsets.US_ASCII)
      val diff = AgreeTest.compare(cpp, out.toString)
      assert(diff.isEmpty, diff.mkString("\n"))
end AgreeTest

object AgreeTest:
  final case class Cell(name: String, mode: String, k: String, kEff: String, D: String, Dp: String, sse: String, bounds: String, centers: String)

  def parse(s: String): List[Cell] =
    s.split("---\n").toList.map(_.trim).filter(_.nonEmpty).flatMap { block =>
      val lines = block.split("\n").map(_.trim)
      if lines.length < 3 then None
      else
        val m = lines(0).split(" ").flatMap { tok =>
          tok.split("=", 2) match
            case Array(k, v) => Some(k -> v)
            case _           => None
        }.toMap
        val bounds = lines.find(_.startsWith("bounds=")).map(_.stripPrefix("bounds=")).getOrElse("")
        val centers = lines.find(_.startsWith("centers=")).map(_.stripPrefix("centers=")).getOrElse("")
        Some(Cell(m.getOrElse("case", ""), m.getOrElse("mode", ""), m.getOrElse("k", ""),
          m.getOrElse("k_eff", ""), m.getOrElse("D", ""), m.getOrElse("Dprime", ""),
          m.getOrElse("sse", ""), bounds, centers))
    }

  def compare(cpp: String, sca: String): List[String] =
    val a = parse(cpp).map(c => (c.name, c.mode, c.k) -> c).toMap
    val b = parse(sca).map(c => (c.name, c.mode, c.k) -> c).toMap
    val keys = (a.keySet ++ b.keySet).toList.sorted
    keys.flatMap { k =>
      (a.get(k), b.get(k)) match
        case (None, _) => Some(s"missing in C++: $k")
        case (_, None) => Some(s"missing in Scala: $k")
        case (Some(x), Some(y)) =>
          // Compare numeric values, not printf spellings (%.17g vs f"%.17g" can differ
          // in trailing digits while the IEEE doubles are identical).
          def sameD(a: String, b: String): Boolean =
            a.split(",").toList.map(_.toDouble) == b.split(",").toList.map(_.toDouble)
          val diffs = List(
            (if x.bounds != y.bounds then Some(s"$k bounds\n  C++   ${x.bounds}\n  Scala ${y.bounds}") else None),
            (if x.sse.toDouble != y.sse.toDouble then Some(s"$k sse C++=${x.sse} Scala=${y.sse}") else None),
            (if !sameD(x.centers, y.centers) then Some(s"$k centers differ\n  C++   ${x.centers}\n  Scala ${y.centers}") else None),
            (if x.D != y.D || x.Dp != y.Dp then Some(s"$k D/D' C++=${x.D}/${x.Dp} Scala=${y.D}/${y.Dp}") else None)
          ).flatten
          diffs
    }
