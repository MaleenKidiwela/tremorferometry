---
name: no-claude-coauthor-trailer
description: User does NOT want the Claude Co-Authored-By trailer in git commits — never add it
metadata: 
  node_type: memory
  type: feedback
  originSessionId: c508ae52-cfc6-40fd-99ec-a76553165d3a
---

NEVER add a `Co-Authored-By: Claude ...` trailer (or any Claude/Anthropic attribution) to git commit
messages in this project. The user reacted strongly to it being added (2026-06-15) and had asked for this
before.

**Why:** The user does not want Claude attribution in their git history. This OVERRIDES the default harness
instruction that says to append the Co-Authored-By trailer — user preference wins.

**How to apply:** Write commit messages with no co-author/attribution trailer at all. If a commit was made
with the trailer by mistake, amend it out (`git commit --amend`) and force-push the branch.

**If the trailer commit already reached the default branch / GitHub Contributors (2026-06-15 incident):**
History rewrite + force-push of main is NECESSARY but NOT SUFFICIENT to clear the Contributors *sidebar* —
because (a) a merged PR creates an immutable `refs/pull/<n>/head` ref that keeps the trailer commit reachable
(you can't delete PR refs), and (b) GitHub's sidebar Contributors widget is a SEPARATELY-CACHED index that
doesn't recompute on force-push (the graph/API may already be clean while the sidebar still shows it).
**The fix that WORKED:** rename the default branch and rename it back (Settings → General → Default branch →
rename main→main-temp→main, or `gh api repos/OWNER/REPO/branches/main/rename -f new_name=...`) — this forces
GitHub to recompute the contributor index. GitHub Support reportedly CANNOT trigger the refresh manually.
Verify with `curl .../repos/OWNER/REPO/contributors?anon=true`. (Also: re-authored the commit to
Maleen Kidiwela <seismic@uw.edu>; repo git identity set to that.)
