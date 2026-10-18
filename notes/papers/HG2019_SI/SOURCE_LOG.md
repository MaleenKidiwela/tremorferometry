# Hutchison & Ghosh (2019) supporting information — source log
Attempted 2026-10-15 ~10:50–11:00 am PT for doi 10.1029/2018JB016138 (Supporting Information S1, Table S1).
- https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB016138 -> HTTP 403, Cloudflare "Just a moment..."
  challenge (`landing.html`); same via the WebFetch tool.
- https://onlinelibrary.wiley.com/doi/full-xml/10.1029/2018JB016138 (Crossref text-mining link) -> 403 challenge.
- https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2018JB016138&file=jgrb53135-sup-0001-...
  and the onlinelibrary.wiley.com equivalent (file name guessed) -> 403 challenge (`s1.bin`, `s2.bin` are the challenge
  pages, not the SI). Wiley TDM API needs a token (no response).
- Crossref metadata (`crossref.json`): no supplement links.
- Unpaywall: only OA copy = ETH Research Collection hdl 20.500.11850/328710 (published version); the handle page is a
  JavaScript app (HTTP 500 to curl), API rate-limited (429/403); sonar.ch record 140124 (`sonar.html`, `sonar.json`) has
  no files attached. figshare search: nothing.
- AGU FM15 abstract (Hutchison & Ghosh, "Toward a complete catalog of VLFEs in Cascadia using a Match Filter Technique",
  studio.m-anage.com/agu/fm15/webprogram/Paper85505.html): early version used "six times the median absolute deviation"
  as threshold and removed time windows with ANSS-catalogued earthquakes; candidates then checked by CMT inversion (low
  VR / inconsistent mechanisms discarded). No station list or counts.
Conclusion: the SI is not accessible without a browser login/challenge. Not used. Ask the user to download it manually
(Table S1 would allow a detection-by-detection match).
- 2026-10-15 11:12 PT: user uploaded Table S1 (csv) and SI S1 (docx) from the publisher; copied here.
