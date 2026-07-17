from __future__ import annotations

from pathlib import Path

from pichia_safe_harbor.annotation_qualification import compare_gff_gene_mirror, summarize_gff, summarize_kegg_list


def test_gff_statistics_preserve_partial_and_functional_counts(tmp_path: Path) -> None:
    gff = tmp_path / "source.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region chr1 1 100\n"
        "chr1\tx\tgene\t1\t10\t.\t+\t.\tID=g1;partial=true;start_range=.,1;end_range=10,.\n"
        "chr1\tx\tmRNA\t1\t10\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tx\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tx\ttRNA\t20\t30\t.\t-\t.\tID=r1\n",
        encoding="utf-8",
    )
    stats = summarize_gff(gff)
    assert stats["gene_count"] == 1
    assert stats["partial_gene_count"] == 1
    assert stats["start_range_gene_count"] == 1
    assert stats["trna_count"] == 1


def test_mirror_and_kegg_identity_are_coordinate_based(tmp_path: Path) -> None:
    left = tmp_path / "left.gff3"
    right = tmp_path / "right.gff3"
    left.write_text("##gff-version 3\n1\tx\tgene\t5\t10\t.\t+\t.\tID=gene-X;Name=X;locus_tag=X\n", encoding="utf-8")
    right.write_text("##gff-version 3\nchr\ty\tgene\t5\t10\t.\t+\t.\tID=gene:X;gene_id=X\n", encoding="utf-8")
    mirror = compare_gff_gene_mirror(left, right)
    assert mirror["shared_gene_id_count"] == 1
    assert mirror["exact_boundary_match_count"] == 0
    mapped = compare_gff_gene_mirror(left, right, {"chr": "1"})
    assert mapped["exact_boundary_match_count"] == 1
    kegg = tmp_path / "kegg.tsv"
    kegg.write_text("ppa:X\tCDS\t1:5..10\tprotein\n", encoding="utf-8")
    assert summarize_kegg_list(kegg)["gene_count"] == 1
