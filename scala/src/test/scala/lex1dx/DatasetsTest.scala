package lex1dx

class DatasetsTest extends munit.FunSuite:

  test("names and unknown"):
    assertEquals(Datasets.names.toList, List("athletics20", "worldathletics-mens-100m-2379"))
    val ex = intercept[IllegalArgumentException] {
      Datasets.load("missing")
    }
    assert(ex.getMessage.contains("unknown dataset"), ex.getMessage)

  test("athletics20"):
    val x = Datasets.load("athletics20")
    assertEquals(x.length, 20)
    assert(x.forall(_.isFinite))
    assert(x.sliding(2).forall(p => p(0) <= p(1)))
    assertEquals(x.head, 9.58)
    assertEquals(x.last, 9.84)
    assertEquals(x.distinct.length, 12)

  test("worldathletics-mens-100m-2379"):
    val x = Datasets.load("worldathletics-mens-100m-2379")
    assertEquals(x.length, 2379)
    assert(x.forall(_.isFinite))
    assert(x.sliding(2).forall(p => p(0) <= p(1)))
    assertEquals(x.head, 9.58)
    assertEquals(x(19), 9.84)
    assert(x.take(20).sameElements(Datasets.load("athletics20")))
    assertEquals(x.distinct.length, 58)

  test("package data files are shipped as resources"):
    def present(name: String): Boolean =
      getClass.getResource(s"/lex1dx/data/$name") != null
    assert(present("athletics20.csv"))
    assert(present("worldathletics-mens-100m-2379.csv"))
    assert(present("README.md"))
end DatasetsTest
