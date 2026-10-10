include somalierpkg/relate

import unittest

proc observation(genotype: int): allele_count =
  case genotype
  of 0: allele_count(nref: 30)
  of 1: allele_count(nref: 15, nalt: 15)
  of 2: allele_count(nalt: 30)
  else: allele_count()

proc sketch(autosomal, x: openArray[int]): counts =
  for genotype in autosomal:
    result.sites.add(observation(genotype))
  for genotype in x:
    result.x_sites.add(observation(genotype))

suite "pairwise genotype counts":
  let workdir = getTempDir() / ("somalier-pairwise-" & $getCurrentProcessId())
  let paths = @[workdir / "a.somalier", workdir / "b.somalier", workdir / "c.somalier"]

  setup:
    createDir(workdir)

  teardown:
    removeDir(workdir)

  test "repeated and reversed lookups preserve autosomal and X counts":
    sketch([0, 1, 2, 0, 1, 2, -1, 0], [0, 1, 2, -1]).write_counts("a", paths[0])
    sketch([2, 1, 2, 1, 0, -1, 1, 0], [2, 1, 2, 0]).write_counts("b", paths[1])
    sketch([-1, -1, -1, -1, -1, -1, -1, -1], [-1, -1, -1, -1]).write_counts("c", paths[2])
    var data = read_extracted(paths, 0.3, 7, false)
    let expected = relation(hets_a: 2, hets_b: 3, hom_alts_a: 2, hom_alts_b: 2,
        ibs0: 1, ibs2: 3, shared_hets: 1, shared_hom_alts: 1, het_ab: 4,
        n: 6, x_ibs0: 1, x_ibs2: 2)
    check data.relatedness(0, 1) == expected
    let missing = data.relatedness(0, 2)
    check missing.n == 0
    check missing.ibs0 == 0
    check missing.ibs2 == 0
    check missing.x_ibs0 == 0
    check missing.x_ibs2 == 0
    check data.relatedness(1, 0) == expected
    check data.relatedness(0, 1) == expected
    check data.relatedness(2, 0) == missing

  test "the combined heterozygous count remains capped at uint16 high":
    let all_het = sketch(newSeqWith(40_000, 1), [])
    all_het.write_counts("a", paths[0])
    all_het.write_counts("b", paths[1])
    var data = read_extracted(paths[0..1], 0.3, 7, false)
    let pair = data.relatedness(0, 1)
    check pair.het_ab == uint16.high
    check pair.n == 40_000
    check pair.shared_hets == 40_000
    check data.relatedness(1, 0) == pair
