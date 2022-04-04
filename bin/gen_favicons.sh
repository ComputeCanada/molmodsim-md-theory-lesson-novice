#!/bin/bash
if [[ $# -lt 1 ]] ; then
    echo "Generate a set of favicons from a high-res template using ImageMagick convert."
    echo "Usage: "
    echo "    $0  path/to/template.png"
    exit 1
elif [[ ! -f $1 ]]; then
    echo "File $1 does not exist."
    exit 1
fi

INPUTFILE="$1"
mkdir -p favicons/

convert $INPUTFILE -resize  57x57    favicons/apple-touch-icon-57x57.png
convert $INPUTFILE -resize  60x60    favicons/apple-touch-icon-60x60.png
convert $INPUTFILE -resize  72x72    favicons/apple-touch-icon-72x72.png
convert $INPUTFILE -resize  76x76    favicons/apple-touch-icon-76x76.png
convert $INPUTFILE -resize 114x114   favicons/apple-touch-icon-114x114.png
convert $INPUTFILE -resize 120x120   favicons/apple-touch-icon-120x120.png
convert $INPUTFILE -resize 144x144   favicons/apple-touch-icon-144x144.png
convert $INPUTFILE -resize 152x152   favicons/apple-touch-icon-152x152.png
convert $INPUTFILE -resize 128x128   favicons/favicon.ico
convert $INPUTFILE -resize  16x16    favicons/favicon-16x16.png
convert $INPUTFILE -resize  32x32    favicons/favicon-32x32.png
convert $INPUTFILE -resize  96x96    favicons/favicon-96x96.png
convert $INPUTFILE -resize 128x128   favicons/favicon-128.png
convert $INPUTFILE -resize 196x196   favicons/favicon-196x196.png
convert $INPUTFILE -resize  70x70    favicons/mstile-70x70.png
convert $INPUTFILE -resize 144x144   favicons/mstile-144x144.png
convert $INPUTFILE -resize 150x150   favicons/mstile-150x150.png
convert $INPUTFILE -resize 310x310   favicons/mstile-310x310.png
convert favicons/mstile-310x310.png -gravity Center -crop 310x150+0+0 +repage  favicons/mstile-310x150.png
