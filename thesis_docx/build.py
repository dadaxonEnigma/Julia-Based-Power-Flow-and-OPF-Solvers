"""
build.py — assemble the final thesis.docx from the Markdown sources using
pandoc and the styled reference.docx (make_reference.py).

  python build.py

Pipeline:
  1. concatenate src/*.md in order
  2. light cleanup of LaTeX→Markdown conversion artifacts (euro symbol,
     equation environments, \\label, broken \\eqref cross-reference links)
  3. pandoc → thesis.docx  (native Word equations, TOC, reference styles)

The Markdown files in src/ are the editable source of truth; this script is
re-runnable so the docx is always reproducible from text.
"""
import os, re, glob
import pypandoc

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
REF = os.path.join(HERE, "reference.docx")
OUT = os.path.join(HERE, "thesis.docx")
BIB = os.path.normpath(os.path.join(HERE, "..", "latex", "literature.bib"))
CSL = os.path.join(HERE, "ieee.csl")   # optional; default style if absent


def cleanup(md: str) -> str:
    # Euro symbol inside math
    md = md.replace(r"\text{\euro}", r"\text{EUR}").replace(r"\euro", "EUR")
    # Drop equation environment wrappers (pandoc keeps display math via $$)
    md = md.replace(r"\begin{equation}", "").replace(r"\end{equation}", "")
    md = md.replace(r"\begin{equation*}", "").replace(r"\end{equation*}", "")
    # Drop \label{...}
    md = re.sub(r"\\label\{[^}]*\}", "", md)
    # Collapse broken eqref/ref links:
    #   [\[eq:x\]](#eq:x){reference-type="eqref" reference="eq:x"} -> (eq.)
    md = re.sub(r"\[\\\[[^\]]*\\\]\]\(#[^)]*\)\{[^}]*\}", "(eq.)", md)
    md = re.sub(r"\[[^\]]*\]\(#[^)]*\)\{reference-type=[^}]*\}", "(ref.)", md)
    # Tidy leftover empty display-math lines
    md = re.sub(r"\$\$\s*\n\s*\n", "$$\n", md)
    return md


def main():
    parts = []
    for path in sorted(glob.glob(os.path.join(SRC, "*.md"))):
        parts.append(cleanup(open(path, encoding="utf-8").read()))
    combined = "\n\n".join(parts)
    # Bibliography rendered under a numbered "References" chapter at the end.
    combined += "\n\n# References\n\n::: {#refs}\n:::\n"

    tmp = os.path.join(HERE, "_combined.md")
    open(tmp, "w", encoding="utf-8").write(combined)

    args = [
        f"--reference-doc={REF}",
        "--toc", "--toc-depth=3",
        "--mathml",                # native Word OMML equations
        "--wrap=preserve",
        "--citeproc",              # resolve [@key] citations + build reference list
        f"--bibliography={BIB}",
    ]
    if os.path.exists(CSL):
        args.append(f"--csl={CSL}")

    pypandoc.convert_file(tmp, "docx", format="markdown", outputfile=OUT, extra_args=args)
    os.remove(tmp)

    # Inject the Innopolis title page (counted but unnumbered) as the first page.
    from title_page import add_title_page
    add_title_page(OUT)

    print("wrote", OUT)


if __name__ == "__main__":
    main()
