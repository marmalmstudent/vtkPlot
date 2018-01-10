all: release
debug:
	sh make_build.sh Debug
release:
	sh make_build.sh Release
clean:
	sh make_clean.sh
