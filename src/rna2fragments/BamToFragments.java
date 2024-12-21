package rna2fragments;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * 
 * Main class
 * 
 */
public class BamToFragments {
	
	public static void printHelp() {
		System.out.println("This software takes aligned 10x or Parse Biosciences RNA-seq data and produces a 10x-style ATAC fragment file.");
		System.out.println("Thus the data can be analyzed using single-cell workflows for ATAC-seq (such as Signac or ArchR).");
		System.out.println("Note that for UMI-deduplication to work, the input file must be position sorted. Only naive deduplication will be performed.");
		System.out.println();
		System.out.println("java -jar rna2fragments.jar gex_possorted_bam.bam  [-m output_match_fragments.tsv] [-n output_intron_fragments.tsv]");
		System.out.println();
		System.out.println("Sort each output:    sort -k 1,1 -k2,2n output_match_fragments.tsv > fragments.sorted.tsv");
		System.out.println("Compress the output: bgzip -@ 8 fragments.sorted.tsv");
		System.out.println("Index the output:    tabix -p vcf fragments.sorted.tsv.gz");
		System.exit(0);
	}
	
	

	/**
	 * 
	 * Entry point
	 * 
	 */
	public static void main(String[] args) {

		//args=new String[]{"/Users/mahogny/Downloads/gex_possorted_bam.bam","gex_fragments.tsv"};

		File foutMatching=null;
		File foutIntron=null;

		
		if(args.length==0) {
			printHelp();
		} else {

			///// Parsing of arguments
			File fin=new File(args[0]);

			int curArg=1;
			while(curArg<args.length) {
				
				if(args[curArg].equals("-m")) {
					foutMatching=new File(args[curArg+1]);
					curArg+=2;
				} else if(args[curArg].equals("-n")) {
					foutIntron=new File(args[curArg+1]);
					curArg+=2;
				} else {
					System.out.println("Parse error on "+args[curArg]);
					System.exit(0);
				}
			}
			
			///// Processing of file
			try {
				PrintWriter pwMatching=null;
				PrintWriter pwIntron=null;
				
				if(foutMatching!=null)
					pwMatching=new PrintWriter(new BufferedWriter(new FileWriter(foutMatching)));
				if(foutIntron!=null)
					pwIntron=new PrintWriter(new BufferedWriter(new FileWriter(foutIntron)));
				
				Processor proc=new Processor();
				proc.read(fin, pwMatching, pwIntron);
				
				if(pwMatching!=null)
					pwMatching.close();
				if(pwIntron!=null)
					pwIntron.close();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
		
	}
	
	
	
	
	
}
