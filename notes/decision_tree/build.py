"""Build the LFE decision-tree page from notes/lfe_decision_tree.json.

  python notes/decision_tree/build.py        -> notes/decision_tree/lfe_decision_tree.html

To update the tree: edit notes/lfe_decision_tree.json (each node: id, parent, kind, status, title,
date, summary, why, lesson, evidence), rebuild, and republish the HTML.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "lfe_decision_tree.json")
KINDS = {"goal", "decision", "catalogue", "attempt", "test", "fix", "result", "future"}
STATUSES = {"success", "partial", "failed", "dead-end", "running", "open", "planned", "superseded"}


def main():
    d = json.load(open(SRC))
    ids = [n["id"] for n in d["nodes"]]
    dup = {i for i in ids if ids.count(i) > 1}
    bad = [n["id"] for n in d["nodes"] if n.get("parent") and n["parent"] not in ids]
    wk = [n["id"] for n in d["nodes"] if n.get("kind") not in KINDS]
    ws = [n["id"] for n in d["nodes"] if n.get("status") not in STATUSES]
    if dup or bad or wk or ws:
        sys.exit(f"invalid tree: duplicate ids {dup}, missing parents {bad}, bad kind {wk}, bad status {ws}")
    tpl = open(os.path.join(HERE, "template.html")).read()
    data = json.dumps(d, ensure_ascii=False).replace("</", "<\\/")
    out = os.path.join(HERE, "lfe_decision_tree.html")
    open(out, "w").write(tpl.replace("/*__DATA__*/", data))
    print(f"wrote {out} ({len(d['nodes'])} nodes)")


if __name__ == "__main__":
    main()
