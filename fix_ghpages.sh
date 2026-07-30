#!/bin/bash
cd .publish-cache/gh-pages && git rm CLAUDE.html && git commit -m "Remove CLAUDE.html" && git push origin HEAD:gh-pages