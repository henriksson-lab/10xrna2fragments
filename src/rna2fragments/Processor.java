package rna2fragments;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class Processor {

	//Statistics of how many records we kept
	int readRecords=0;
	int keptRecords=0;
	String currentSource=null;
	int unalignedRecord=0;
	int skippedDup=0;
	int skippedWrongBC=0;
	int keptBlocks=0;

	
	public void read(File fBAM, PrintWriter pwMatching, PrintWriter pwIntron) throws IOException {
		

		//Open the BAM file and get to work
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		
		

		
		//This is to keep track of duplicates.
		//Approximate, as the same cDNA can be fragmented multiple times in the library prep
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
				
		//Loop through all SAM records
		for (final SAMRecord samRecord : reader) {

			//Get UMI and BC for this read
			String bcCellCurrentUMI=getUMI(samRecord);
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

					//Only consider aligned reads i.e. those having a contig reference
					if(samRecord.getContig()!=null) {
						
						
						//Print matching parts of CIGAR
						if(pwMatching!=null) {

							//Iterate over blocks that have been directly aligned to the reference sequence - from CIGAR string
							//same as: SAMUtils.getAlignmentBlocks(getCigar(), getAlignmentStart(), "read cigar");
							List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
							for(int curAB=0;curAB<listBlocks.size();curAB++) {
								AlignmentBlock ab=listBlocks.get(curAB);
								
								String blockSource=samRecord.getContig();
								int blockFrom=ab.getReferenceStart();
								int blockTo=ab.getReferenceStart()+ab.getLength();

								printBlock(pwMatching, bcCellCurrentCellBarcode, 
										blockSource,blockFrom,blockTo);
								
								keptBlocks++;							
							}

						}
						

						//Print intron-like parts of CIGAR ("N")
						if(pwIntron!=null) {
							List<IntronBlock> listBlocks=getIntronBlocks(samRecord.getCigar(), samRecord.getAlignmentStart(), "read cigar");;
							for(int curAB=0;curAB<listBlocks.size();curAB++) {
								IntronBlock ab=listBlocks.get(curAB);
								
								String blockSource=samRecord.getContig();
								int blockFrom=ab.getReferenceStart();
								int blockTo=ab.getReferenceStart()+ab.getLength();

								printBlock(pwIntron, bcCellCurrentCellBarcode, 
										blockSource,blockFrom,blockTo);
								
								keptBlocks++;							
							}

						}
						
						keptRecords++;
					} else {
						//Unaligned read
						unalignedRecord++;
					}
					
				}
			} else {
				skippedWrongBC++;
			}
			
			//Print current status after this record has been processed
			readRecords++;
			if(readRecords%1000000 == 0)
				printStatistics();

		}
	

	//Print final statistics
		printStatistics();

	reader.close();
	}
	
	
	public void printStatistics() {
		System.out.println(
				"Kept rec/Read rec: "+keptRecords+"/"+readRecords+
				"  @sequence: "+currentSource+
				"  BadBC: "+skippedWrongBC+
				"  Unaligned: "+unalignedRecord+
				"  Blocks written: "+keptBlocks+
				"  SkipDup: "+skippedDup);
	}
	
	public static void printBlock(PrintWriter pw, String bcCellCurrentCellBarcode, String blockSource, int blockFrom, int blockTo) {
		pw.println(blockSource+"\t"+blockFrom+"\t"+blockTo+"\t"+bcCellCurrentCellBarcode+"\t1");
	}
	

	public static String getUMI(SAMRecord samRecord) {
		return (String)samRecord.getAttribute("UB");

		
		//TODO: flag to support PB UMIs; need to be parsed differently

		//this is the name of a Parse Bioscience's read:
		//01_01_35__T__1_1_35__CATTCCTA_AACGTGAT_CATACCAA__TTCGCTCCAA__241218IC__OH_@LH00202:171:22JJK7LT4:3:2208:30933:24827
		//UMI is thus part of the name, but it is not a tag!
		
		//TODO: option to store introns!
		
	}


   /**
    * Given a Cigar, Returns blocks of the sequence that seem to be intronic (N, skips against reference)
    * 
    * modified from htslib, MIT license
    *
    * @param cigar          The cigar containing the alignment information
    * @param alignmentStart The start (1-based) of the alignment
    * @param cigarTypeName  The type of cigar passed - for error logging.
    * @return List of alignment blocks
    */
   public static List<IntronBlock> getIntronBlocks(final Cigar cigar, final int alignmentStart, final String cigarTypeName) {
       if (cigar == null) return Collections.emptyList();

       final List<IntronBlock> alignmentBlocks = new ArrayList<>();
       int readBase = 1;
       int refBase = alignmentStart;

       for (final CigarElement e : cigar.getCigarElements()) {
           switch (e.getOperator()) {
               case H:
                   break; // ignore hard clips
               case P:
                   break; // ignore pads
               case S:
                   readBase += e.getLength();
                   break; // soft clip read bases
               case N:
                   final int lengthN = e.getLength();
                   refBase += lengthN;
                   alignmentBlocks.add(new IntronBlock(readBase, refBase, lengthN));
                   break;  // reference skip
               case D:
                   refBase += e.getLength();
                   break;
               case I:
                   readBase += e.getLength();
                   break;
               case M:
               case EQ:
               case X:
                   final int length = e.getLength();
                   readBase += length;
                   refBase += length;
                   break;
               default:
                   throw new IllegalStateException("Case statement didn't deal with " + cigarTypeName + " op: " + e.getOperator() + "in CIGAR: " + cigar);
           }
       }
       return Collections.unmodifiableList(alignmentBlocks);
   }

	
	
	
}