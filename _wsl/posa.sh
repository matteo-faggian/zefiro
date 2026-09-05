#!/bin/bash
# Porta il repo ricostruito dentro la cartella di lavoro, e ci mette dentro
# anche le ultime modifiche (gitattributes, gitignore, script nuovi).
set -e
R=/root/zefiro_git
M=/mnt/i/AA_ENGINE
[ -d "$R/.git" ] || { echo "il repo ricostruito non c'e'"; exit 2; }
[ -d "$M/.git" ] && { echo "ATTENZIONE: I:\\AA_ENGINE ha gia' un .git, non tocco niente"; exit 3; }
cp -r "$R/.git" "$M/.git"
cd "$M"
git config core.autocrlf false
git config user.name "mat"
git config user.email "gptmatandrew@gmail.com"
echo "=== che cosa risulta cambiato rispetto al commit finale ==="
git status --short | head -20
echo
git add -A
git -c user.name=mat -c user.email=gptmatandrew@gmail.com commit -q -m "$(cat <<'MSG'
Regole di fine-riga e quarantena esclusa dal versionamento

`.gitattributes`: nel repo si sta in LF sempre. Il master e' su un disco
Windows e il calcolo gira in WSL: senza questa regola basta un editor che
salva CRLF perche' il file dopo risulti riscritto da capo nel diff, e perche'
bash si lamenti di caratteri invisibili a fine riga. Oggi sono 241 file su 241
in LF: la regola serve a tenerceli.

`.gitignore`: esclusa la cartella (eliminare), che e' quarantena e non
progetto. Dentro pero' ci sono ancora i quattordici zip datati da cui e' stata
ricostruita la storia di questo repo: finche' il push non e' verificato, quelli
sono l'unica copia e non vanno svuotati.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PqPTkjhc82NmVTJozaYyCp
MSG
)"
echo "=== stato finale ==="
git log --oneline | head -3
echo "..."
git log --oneline | tail -2
echo
echo "commit totali: $(git rev-list --count HEAD)"
echo "file versionati: $(git ls-tree -r HEAD --name-only | wc -l)"
git status --short | head -5
echo "(niente sopra = albero pulito)"
