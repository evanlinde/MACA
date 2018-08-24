#!/bin/bash
#
# Convert the Jupyter notebook (master document) into HTML for
# web display and print-friendly PDFs.
# 

file_basename="$1"

declare -A titles
titles=(
  [MACA_to_Envision_SWAT]="\"Converting MACA Data for Envision and SWAT Models\""
  [MACA_to_SWAT]="\"Converting MACA Data for SWAT Models\""
  [AWR]="\"Converting MACA Data for Envision and SWAT Models - AWR\""
)
declare -A modify_dates
modify_dates=(
  [MACA_to_Envision_SWAT]="2017-10-19"
  [MACA_to_SWAT]="2017-11-21"
  [AWR]="2018-08-24"
)
declare -A pdf
pdf=(
  [MACA_to_Envision_SWAT]=1
  [MACA_to_SWAT]=1
  [AWR]=0
)

#title="\"Converting MACA Data for Envision and SWAT Models\""
title=${titles[${file_basename}]}
author="\"Evan Linde\""
institution="\"Oklahoma State University\""
#today="\"$(date +'%B %-d, %Y')\""
[[ ${modify_dates[${file_basename}]+1} ]] && today="\"$(date --date="${modify_dates[${file_basename}]}" +'%B %-d, %Y')\"" || today="\"$(date +'%B %-d, %Y')\""
#file_basename="MACA_to_Envision_SWAT"

# Convert the notebook to markdown
# The markdown file will be ${file_basename}.md
jupyter nbconvert --to markdown ${file_basename}.ipynb

# Fix image references in the markdown file. 
# The image references in the notebook are all formatted like
#    ![relative/path/to/image.png](attachment:image.png)
# so we're replacing what's in the parentheses with what's
# in the square brackets.
sed -i 's/^\!\[\(.*\)\](attachment:.*)\ *$/\!\[\1\]\(\1\)/' ${file_basename}.md

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

# Then add the body from the exported markdown document
# fixing image references as we go. 
# The image references in the notebook are all formatted like
#    ![relative/path/to/image.png](attachment:image.png)
# so we're replacing what's in the parentheses with what's
# in the square brackets.
sed -n '/^\#\ Introduction/,${s/^\!\[\(.*\)\](attachment:.*)\ *$/\!\[\1\]\(\1\)/;p}' ${file_basename}.md >> doc/${file_basename}.md

# Get rid of the directly exported markdown since we're done with it
rm ${file_basename}.md

# Make HTML page
pandoc -f markdown_github+tex_math_dollars+yaml_metadata_block doc/${file_basename}.md --css doc/pandoc.css --toc -s --mathjax -o ${file_basename}.html
# Load mathjax over https instead of http
sed -i 's,http://cdn.mathjax.org,https://cdn.mathjax.org,g' ${file_basename}.html


# Only make PDFs for some documents
if [[ ${pdf[${file_basename}]} -eq 1 ]]; then
  # Make a PDF with default margins
  pandoc -f markdown_github+tex_math_dollars+yaml_metadata_block --toc --listings -H doc/listings-setup.tex doc/${file_basename}.md -o doc/${file_basename}.pdf

  # Make a PDF with smaller margins
  pandoc -f markdown_github+tex_math_dollars+yaml_metadata_block --toc --listings -H doc/listings-setup.tex -V geometry:"left=2.5cm, top=2cm, right=2.5cm, bottom=2cm" doc/${file_basename}.md -o doc/${file_basename}_small_margins.pdf
fi
