# Observatory links, 2026-10-09

Pulled the local RiboCrypt master checkout to `3c72c86` from
`https://github.com/m-swirski/RiboCrypt`. HTTPS was used because this execution
environment has no SSH executable. No backend source edits were made.

The current API still uses `#Observatory?obs_state=...` and a version-1 compressed
state. Its restore logic validates controls and requires the saved cohorts for
browser autostart. The dominant-state app and review-link generator duplicated
that codec and appended an empty `All merged` group as a previous restore guard.
That group can expand to the complete collection, and expanded membership does
not match the explicitly saved empty cohort. It also introduces an unwanted
whole-collection comparison track.

The app and link generator now share `scripts/dominant_observatory_urls.R`, which
loads the local development RiboCrypt with `devtools::load_all()` and delegates
state construction, validation, parsing, and encoding to its URL functions.
Matched links contain their actual case/control groups, with a populated active
group. Partial-pooling groups are retained. Empty All merged guard groups are
removed; genuinely populated All merged groups are preserved. Links without a
matched cohort open the transcript with autostart disabled for manual selection.

Older links in review, label, and workbench tables are normalized when the app
reads them, preserving the hostname/application path, transcript, controls,
cohort membership, and manual labels. The app also finds the metadata CSV under
`~/livemount/Bio_data/NGS_pipeline/` when the old S-drive path is unavailable.
The generated review CSVs, HTML/Markdown cards, and static scientific-note docs
were refreshed. Pipeline dependencies now include the backend URL and restore
helper files so future API changes invalidate generated review links.

Validation:
- 228 package assertions passed with no test warnings or failures.
- App startup completed and all 946 distinct card/workbench URLs parsed through
  the pulled API; none retained an empty All merged guard.
- Focused tests cover the old guard/expanded-cohort readiness mismatch, corrected
  autostart readiness, exact cohort preservation, partial/populated selections,
  missing-cohort handling, URL-path preservation, idempotent upgrades, and manual
  review-label preservation.
- `git diff --check` passed. Logs are under
  `results/observatory_url_fix_2026-10-09/`.

The live public Observatory was not driven through a browser in this check.
The verification uses its current local parser and restore/readiness functions.
An already running dominant-state app must reload its application code to pick
up these changes; restart that app process if refreshing alone retains old links.
