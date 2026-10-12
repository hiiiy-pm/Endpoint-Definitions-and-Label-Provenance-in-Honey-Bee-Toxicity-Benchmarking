"""Check actual transformed font sizes at the manuscript's figure scale.

The smaller PDF is a QA facsimile of embedding, not a replacement figure.
The source table and archived comparison copy must remain byte-identical.
"""
from pathlib import Path
import hashlib
import json
import pymupdf as fitz

QA = Path(__file__).resolve().parent
ROOT = QA.parents[2]
SCALE = .712


def spans(page):
    return [span for block in page.get_text("dict")["blocks"] if "lines" in block
            for line in block["lines"] for span in line["spans"] if span["text"].strip()]


def main():
    source_path = ROOT / "figures/Fig7_label_crossover.pdf"
    source_csv = ROOT / "source_data/figure_label_crossover.csv"
    old_csv = QA / "before_readability_repair/figure_label_crossover.csv"
    assert source_csv.read_bytes() == old_csv.read_bytes(), "Source table changed during presentation work"
    source = fitz.open(source_path)
    page = source[0]
    original = spans(page)
    assert min(s["size"] for s in original) >= 10 - 1e-5
    smaller = fitz.open()
    target = smaller.new_page(width=page.rect.width * SCALE, height=page.rect.height * SCALE)
    target.show_pdf_page(target.rect, source, 0)
    output = QA / "Fig7_manuscript_scale_0712.pdf"
    smaller.save(output, deflate=True)
    smaller.close()
    smaller = fitz.open(output)
    transformed = spans(smaller[0])
    assert len(transformed) == len(original)
    minimum = min(s["size"] for s in transformed)
    assert minimum >= 7.0, minimum
    for span in transformed:
        assert smaller[0].rect.contains(fitz.Rect(span["bbox"])), span["text"]
    smaller[0].get_pixmap(dpi=144, alpha=False).save(QA / "Fig7_manuscript_scale_0712_preview.png")
    report = {"Status": "passed", "embedding_scale": SCALE,
              "source_inches": [page.rect.width / 72, page.rect.height / 72],
              "embedded_inches": [smaller[0].rect.width / 72, smaller[0].rect.height / 72],
              "source_minimum_pt": min(s["size"] for s in original), "embedded_minimum_pt": minimum,
              "embedded_maximum_pt": max(s["size"] for s in transformed),
              "text_spans_checked": len(transformed), "all_text_inside_page": True,
              "source_csv_unchanged": True,
              "source_csv_sha256": hashlib.sha256(source_csv.read_bytes()).hexdigest(),
              "source_pdf_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
              "limitation": "Checks the declared 0.712 embedding scale; recheck the compiled manuscript if its placement scale changes."}
    (QA / "MANUSCRIPT_SCALE_QA.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
