all: favicon-114x114.png favicon-72x72.png favicon-54x54.png favicon.ico

clean:
	rm -f favicon-114x114.png favicon-72x72.png favicon-54x54.png favicon.ico

favicon-114x114.png: icon.png
	convert -geometry 114x114 icon.png favicon-114x114.png

favicon-72x72.png: icon.png
	convert -geometry 72x72 icon.png favicon-72x72.png

favicon-54x54.png: icon.png
	convert -geometry 54x54 icon.png favicon-54x54.png

favicon.ico: favicon.png
	convert favicon.png -bordercolor white -border 0 \
          \( -clone 0 -resize 16x16 \) \
          \( -clone 0 -resize 32x32 \) \
          \( -clone 0 -resize 48x48 \) \
          \( -clone 0 -resize 64x64 \) \
          -delete 0 -alpha off -colors 256 favicon.ico
