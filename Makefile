all:
	cd src; javac -cp .:../htsjdk-2.23.0-3-g657b0a6-SNAPSHOT.jar rna2fragments/BamToFragments.java 
	cd src; jar cfm ../rna2fragments.jar Manifest.txt ./rna2fragments/*.class

gitaddall:
	git add */*.java

