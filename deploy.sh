#!/bin/bash
# deploy.sh — deploy LociTracker_site to GitHub Pages (gh-pages branch)

# -----------------------------
# CONFIG
# -----------------------------
SITE_FOLDER="LociTracker_site"
BRANCH="gh-pages"
REMOTE="origin"

# -----------------------------
# 1. Make sure we are on main
# -----------------------------
git checkout main || { echo "Failed to checkout main"; exit 1; }

# -----------------------------
# 2. Commit any changes to main
# -----------------------------
git add "$SITE_FOLDER"
git commit -m "Update ShinyLive site source" 2>/dev/null

# -----------------------------
# 3. Switch/create gh-pages branch
# -----------------------------
git checkout $BRANCH 2>/dev/null || git checkout -b $BRANCH

# -----------------------------
# 4. Delete everything in gh-pages branch
# -----------------------------
git rm -rf . > /dev/null 2>&1

# -----------------------------
# 5. Copy contents of site folder to root
# -----------------------------
cp -R ../$SITE_FOLDER/* . 2>/dev/null || cp -R $SITE_FOLDER/* .

# -----------------------------
# 6. Add all files and commit
# -----------------------------
git add .
git commit -m "Deploy ShinyLive site" 2>/dev/null

# -----------------------------
# 7. Force push to remote
# -----------------------------
git push $REMOTE $BRANCH --force

# -----------------------------
# 8. Switch back to main
# -----------------------------
git checkout main

echo "Deployment complete! Your site should be live at:"
echo "https://<username>.github.io/<repo>/"
