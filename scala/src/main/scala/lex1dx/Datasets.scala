package lex1dx

import scala.io.Source

/** Bundled example datasets.

  * Source: World Athletics 100m all-time lists, retrieved January 2024.
  *
  * Sample Data Attribution: The sample 100m timing datasets (top 20 performances and 2379 all time bests) contain historical athletic performance records compiled from publicly available statistics published by World Athletics (https://worldathletics.org/). These datasets are provided strictly as static benchmark/example inputs for testing and demonstrating the algorithms in this package. This project is not affiliated with, endorsed by, or sponsored by World Athletics.
  */
object Datasets:
  def names: IArray[String] = IArray("athletics20", "worldathletics-mens-100m-2379")

  def load(name: String): Array[Double] =
    val file = name match
      case "athletics20"                       => "athletics20.csv"
      case "worldathletics-mens-100m-2379"     => "worldathletics-mens-100m-2379.csv"
      case _ =>
        throw new IllegalArgumentException(
          s"unknown dataset '$name'; known names: ${names.mkString(", ")}"
        )
    val path = s"/lex1dx/data/$file"
    val in = Option(getClass.getResourceAsStream(path)).getOrElse(
      throw new IllegalStateException(s"missing dataset resource $path")
    )
    val src = Source.fromInputStream(in, "UTF-8")
    try
      val rows = src.getLines().zipWithIndex.flatMap { (line, i) =>
        val t = line.strip
        if t.isEmpty || (i == 0 && t.toLowerCase == "time") then None
        else Some(t.toDouble)
      }.toArray
      if rows.isEmpty then throw new IllegalArgumentException(s"dataset '$name' is empty")
      rows
    finally src.close()
end Datasets
