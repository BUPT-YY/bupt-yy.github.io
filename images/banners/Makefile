all: teaching.jpg software.jpg home.jpg

teaching.jpg: Chalkandboard.jpg Makefile
	convert $< -gravity Northwest -crop 75% +repage -resize 1280x1280 -fill black -colorize 50% $@

home.jpg: sarah-and-jim.jpg Makefile
	convert $< -gravity North -crop 80% +repage -resize 1280x1280 -fill black -colorize 40% $@

software.jpg: tasha-and-calculus.jpg Makefile
	convert $< -gravity South -crop 75% +repage -resize 1280x1280 -gaussian-blur 5x5 -fill black -colorize 30% $@
