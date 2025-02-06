<!--
SPDX-PackageName: "detray, a part of the ACTS project"
SPDX-FileCopyrightText: 2021 CERN
SPDX-License-Identifier: MPL-2.0
-->

# Google Benchmark Build Instructions

This subdirectory holds instructions for building
[benchmark](https://github.com/google/benchmark) as part of this project.
This is meant to come in handy for building the project's benchmarks in
environments which do not provide Google Branchmark themselves.

Note that since Google Benchmark is only needed for the tests of this project,
which are not installed together with the project, Google Benchmark is not
installed together with the project either.
