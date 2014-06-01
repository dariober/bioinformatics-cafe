package bisReadBias;

import net.sf.samtools.*;
import net.sf.picard.reference.*;

import java.io.*;

public class Main {

	/**
	 * Return the read index where the alignment starts after having skipped
	 * soft and hard clipped bases.
	 * E.g. if the read aligns with cigar 70M, then
	 *     getReadStartEndIndexes(rec) >>> 0, 70
	 * If cigar 65MS5:
	 * 	   getReadStartEndIndexes(rec) >>> 0, 65
	 * 
	 * The returned indexes can be used to iterate through the reference
	 * bases correpsonding to each aligned read base.   
	 * 
	 * @param rec
	 * @return Array of length 2: {Index start, Index end}
	 */
	private static int[] getReadStartEndIndexes(SAMRecord rec){

		int alignedReadStartIdx= rec.getAlignmentStart() - rec.getUnclippedStart();
		int alignedReadEndIdx= rec.getReadLength() - (rec.getUnclippedEnd() - rec.getAlignmentEnd());

		int[] startEnd= {alignedReadStartIdx, alignedReadEndIdx};
		return(startEnd);
	}
	
	public static void main(String[] args) {

		String fastaRef= "test_data/chr7.fa";
		String bam= "test_data/ds051.actb.clip.sam";
		
		IndexedFastaSequenceFile fastaFile= null;
		try {
			fastaFile= new IndexedFastaSequenceFile(new File(fastaRef));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		SAMFileReader samfile= new SAMFileReader(new File(bam));

		for (SAMRecord rec : samfile){
			
			int[] startEnd= getReadStartEndIndexes(rec);			
			int alignedReadStartIdx= startEnd[0];
			int alignedReadEndIdx= startEnd[1];
			
			byte readbases[]= rec.getReadBases();
			for (int readidx= alignedReadStartIdx; readidx < alignedReadEndIdx; readidx++){
				
				byte b= readbases[readidx];
				
				StringBuilder sb= new StringBuilder();
				int refidx= rec.getReferencePositionAtReadPosition(readidx);
				byte[] refb= fastaFile.getSubsequenceAt(rec.getReferenceName(), refidx+1, refidx+1).getBases();
			
				String ref = new String(refb);

				sb.append( (readidx+1) + " ");
				sb.append( (char)b + " " );
				sb.append(ref + " ");				
				System.out.println(sb.toString());
			}
			break;
		} // End loop through sam records
		
		// Close 
		samfile.close();
		try {
			fastaFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}