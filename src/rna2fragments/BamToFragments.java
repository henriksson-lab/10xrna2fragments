package rna2fragments;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;


import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamToFragments {
	
	
	public static void read(File fBAM, PrintWriter pw) throws IOException {

		

		//Open the BAM file and get to work
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		
		
		//Statistics of how many records we kept
		int readRecords=0;
		int keptRecords=0;
		String currentSource=null;
		int unalignedRecord=0;
		int skippedDup=0;
		int skippedWrongBC=0;

		
		//This is to keep track of duplicates.
		//Approximate, as the same cDNA can be fragmented multiple times in the library prep
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
				
		//Loop through all SAM records
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%1000000 == 0){
				
				System.out.println(
						"Kept/Read: "+keptRecords+"/"+readRecords+
						"  @sequence: "+currentSource+
						"  BadBC: "+skippedWrongBC+
						"  Unaligned: "+unalignedRecord+
						"  SkipDup: "+skippedDup);

				
			}

			//Get UMI and BC for this read
			String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
			String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
				
			//If the read has no BC then ignore it
			if(bcCellCurrentCellBarcode!=null) {
				//Check if duplicate read, if UMI present; ATAC, dedup by coordinate?
				if(bcCellCurrentUMI!=null && bcCellCurrentUMI.equals(bcCellPreviousU) && bcCellCurrentCellBarcode.equals(bcCellPreviousC)) {
					//Do nothing, just ignore - this is a duplicate read
					skippedDup++;
				} else {
					//Remember for later
					bcCellPreviousU=bcCellCurrentUMI;
					bcCellPreviousC=bcCellCurrentCellBarcode;

					//Store position, for displaying progress
					currentSource=samRecord.getContig();

					
					if(samRecord.getContig()!=null) {
						///////////////// For each block
						
						//A read may have been split into multiple blocks. 
						//Count these separately. Naive assumption that these are split over introns... is this correct?
						List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
						for(int curAB=0;curAB<listBlocks.size();curAB++) {
							AlignmentBlock ab=listBlocks.get(curAB);
							
							String blockSource=samRecord.getContig();
							int blockFrom=ab.getReferenceStart();
							int blockTo=ab.getReferenceStart()+ab.getLength();

							
							pw.println(blockSource+"\t"+blockFrom+"\t"+blockTo+"\t"+bcCellCurrentCellBarcode+"\t1");
							keptRecords++;
							
						}
						
					} else {
						//Unaligned read
						unalignedRecord++;
					}
					
				}
			} else {
				skippedWrongBC++;
			}
		}
	

		
	System.out.println(
			"Kept/Read: "+keptRecords+"/"+readRecords+
			"  @sequence: "+currentSource+
			"  BadBC: "+skippedWrongBC+
			"  Unaligned: "+unalignedRecord+
			"  SkipDup: "+skippedDup);

	reader.close();
		
	}
	
	
	
	public static void printHelp() {
		System.out.println("This software takes 10x rnaseq and produces a 10x-like fragment file.");
		System.out.println("Thus they can be analyzed using single-cell workflows (such as signac or archr)");
		System.out.println();
		System.out.println("java -jar 10xrna2fragments.jar gex_possorted_bam.bam gex_fragments.tsv");
		System.out.println();
		System.out.println("Sort it: sort -k 1,1 -k2,2n gex_fragments.tsv > gex_fragments.sorted.tsv");
		System.out.println("Compress it: bgzip -@ 8 gex_fragments.sorted.tsv");
		System.out.println("Index it: tabix -p vcf gex_fragments.sorted.tsv.gz");
		System.exit(0);
	}
	
	
	public static void main(String[] args) {

		//args=new String[]{"/Users/mahogny/Downloads/gex_possorted_bam.bam","gex_fragments.tsv"};
		
		if(args.length==0) {
			printHelp();
		}
		
		if(args.length==2) {
			
			File fin=new File(args[0]);
			File fout=new File(args[1]);
			
			try {
				PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(fout)));
				
				read(fin, pw);
				
				pw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		} else {
			System.out.println("Parse error on "+args[0]);
			System.exit(0);
		}
		
		
		
	}


}
