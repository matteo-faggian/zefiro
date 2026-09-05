#!/bin/bash
# Completa la posa del repo. Il .git e' gia' copiato; mancava solo dire a git
# che una cartella su un disco Windows, vista da root in WSL, e' sua.
# (E' la protezione contro i repo altrui in cartelle condivise: qui il
#  proprietario "diverso" e' solo l'artefatto del mount DrvFs.)
set -e
M=/mnt/i/AA_ENGINE
git config --global --add safe.directory "$M"
cd "$M"
git config core.autocrlf false
git config user.name "mat"
git config user.email "gptmatandrew@gmail.com"
echo "=== cambiato rispetto al commit finale ==="
git status --short | head -20
echo
git add -A
git commit -q -F - <<'MSG'
Regole di fine-riga e quarantena esclusa dal versionamento

`.gitattributes`: nel repo si sta in LF sempre. Il master e' su un disco
Windows e il calcolo gira in WSL: senza questa regola basta un editor che
salva CRLF perche' il file dopo risulti riscritto da capo nel diff, e perche'
bash si lamenti di caratteri invisibili a fine riga. Oggi sono 241 file su 241
in LF: la regola serve a tenerceli.

`.gitignore`: esclusa la cartella (eliminare), che e' quarantena e non
progetto. Dentro pero' ci sono ancora i quattordici zip datati da cui e' stata
ricostruita la storia di questo repo: finche' il push non e' verificato,
quelli sono l'unica copia e non vanno svuotati.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PqPTkjhc82NmVTJozaYyCp
MSG
echo "=== stato finale ==="
git log --oneline | head -3
echo "   ..."
git log --oneline | tail -2
echo
echo "commit totali : $(git rev-list --count HEAD)"
echo "file versionati: $(git ls-tree -r HEAD --name-only | wc -l)"
echo "peso          : $(du -sh .git | cut -f1)"
echo "--- albero pulito? ---"
git status --short | head -5
echo "(niente sopra = pulito)"
