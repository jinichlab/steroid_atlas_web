#!/bin/bash
# Wait for score audit to finish, then classify + run master-table assembly.
# Fully self-contained — launch under nohup + disown and it survives any laptop close.
set -e
cd /home/adsiordia/steroid_atlas
export LD_LIBRARY_PATH=/home/adsiordia/miniconda3/lib

PY=/home/adsiordia/miniconda3/bin/python3
LOG=/tmp/final_master_chain.log

echo "=== $(date) === waiting for score audit to finish..." >> $LOG
# Wait until the score audit process is gone
while pgrep -f "15_annotation_score_audit" >/dev/null; do
    sleep 5
done
echo "=== $(date) === score audit finished, running classification..." >> $LOG

# Re-run the score audit's classification step (uses full cache, skips retry)
$PY <<'PY_INNER' >> $LOG 2>&1
import json, pandas as pd
cache = json.load(open('analysis/score_audit_cache.json'))
p = pd.read_csv('data/proteins.csv', low_memory=False)
rows = []
for _, r in p.iterrows():
    acc = r['accession']
    ann = cache.get(acc, {})
    etype = ann.get('entry_type','') or ''
    score = float(ann.get('annotation_score', 0) or 0)
    reason = ann.get('inactive_reason','') or ''
    if etype.startswith('Inactive') or 'DELETED' in reason.upper() or 'Not part' in reason:
        decision = 'DROP_INACTIVE'
    elif not ann:
        decision = 'DROP_UNFETCHED'
    elif score < 5.0:
        decision = 'DROP_LOWSCORE'
    elif score == 5.0:
        decision = 'KEEP'
    else:
        decision = 'REVIEW'
    rows.append({'accession': acc, 'entry_type': etype,
                 'annotation_score': score, 'inactive_reason': reason,
                 'decision': decision})
audit = pd.DataFrame(rows)
audit.to_csv('analysis/score_audit.tsv', sep='\t', index=False)
print(f'classified {len(audit):,} entries')
print(audit['decision'].value_counts().to_string())
PY_INNER

echo "=== $(date) === running master table assembly..." >> $LOG
$PY analysis/16_fetch_pubmed_and_build_master.py >> $LOG 2>&1

echo "=== $(date) === DONE. Master table written to data/proteins.csv" >> $LOG
