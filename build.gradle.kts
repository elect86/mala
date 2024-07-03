plugins {
    embeddedKotlin("jvm")
}

repositories { mavenCentral() }

dependencies {
    implementation("ai.djl:api")
    implementation("ai.djl.pytorch:pytorch-engine")
    implementation("ai.djl:djl.kt")
    implementation("org.ejml:ejml-all:0.43")

    implementation("org.jetbrains.kotlinx:multik-core:0.2.3")
    implementation("org.jetbrains.kotlinx:multik-default:0.2.3")

    implementation(kotlin("test"))
}

tasks {
    test { useJUnitPlatform() }
}