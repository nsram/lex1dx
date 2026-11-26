ThisBuild / organization := "io.github.nsram"
ThisBuild / version := "0.1.3"
ThisBuild / scalaVersion := "3.3.4"
ThisBuild / licenses := Seq("MIT" -> url("https://opensource.org/licenses/MIT"))

lazy val root = (project in file("."))
  .settings(
    name := "lex1dx",
    description := "lex1dx: fast 1-D k-means with a 0.01% REE tolerance",
    scalacOptions ++= Seq("-deprecation", "-feature", "-Wunused:all"),
    libraryDependencies += "org.scalameta" %% "munit" % "1.0.0" % Test,
    Test / fork := true,
    Test / javaOptions += "-Xmx2g",
    Test / baseDirectory := (ThisBuild / baseDirectory).value
  )
