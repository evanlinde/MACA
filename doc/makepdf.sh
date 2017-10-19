#!/bin/bash

# This script is intended to be called from outside the doc directory
# i.e. where the notebook file lives. For example:
#
#    bash doc/makepdf.sh
# 

title="\"Converting MACA Data for Envision and SWAT Models\""
author="\"Evan Linde\""
institution="\"Oklahoma State University\""
today="\"$(date +'%B %-d, %Y')\""
file_basename="MACA_to_Envision_SWAT"

# Convert the notebook to markdown
# The markdown file will be ${file_basename}.md
jupyter nbconvert --to markdown ${file_basename}.ipynb

# Build the markdown file we'll feed to pandoc to make the PDF
# Start with the header
cat << EOF > doc/${file_basename}.md
---
title: ${title}
author: ${author}
institution: ${institution}
date: ${today}
---


EOF

# Then add the body from the exported markdown document while
# fixing image references as they're encountered.
# The image references in the notebook are all formatted like
#    ![/full/path/to/image.png](attachment:image.png)
# so we're replacing what's in the parentheses with what's
# in the square brackets.
sed -n '/^\#\ Introduction/,${s/^\!\[\(.*\)\](attachment:.*)\ *$/\!\[\1\]\(\1\)/;p}' ${file_basename}.md >> doc/${file_basename}.md

rm ${file_basename}.md

# Make a PDF with default margins
pandoc -f markdown_github+tex_math_dollars+yaml_metadata_block --toc --listings -H doc/listings-setup.tex doc/${file_basename}.md -o doc/${file_basename}.pdf

# Make a PDF with smaller margins
pandoc -f markdown_github+tex_math_dollars+yaml_metadata_block --toc --listings -H doc/listings-setup.tex -V geometry:"left=2.5cm, top=2cm, right=2.5cm, bottom=2cm" doc/${file_basename}.md -o doc/${file_basename}_small_margins.pdf

