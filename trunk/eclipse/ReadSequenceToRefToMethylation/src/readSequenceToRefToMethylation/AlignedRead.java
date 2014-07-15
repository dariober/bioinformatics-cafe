package readSequenceToRefToMethylation;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.biojava3.core.sequence.DNASequence;





// import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * @author berald01
 * 
 * TODO: Define a CompoundSet where soft clipped is marked by a special character.
 * 
 * MEMO on cigar oprators: 
 * -----------------------
 * 		System.out.println("Read base consumers");
		System.out.println(CigarOperator.M.consumesReadBases());  //True
		System.out.println(CigarOperator.I.consumesReadBases());
		System.out.println(CigarOperator.EQ.consumesReadBases());
		System.out.println(CigarOperator.S.consumesReadBases());
		System.out.println(CigarOperator.X.consumesReadBases());
		
		System.out.println(CigarOperator.D.consumesReadBases()); //False
		System.out.println(CigarOperator.N.consumesReadBases());
		System.out.println(CigarOperator.P.consumesReadBases());
		System.out.println(CigarOperator.H.consumesReadBases());
		
		System.out.println("Reference base consumers");
		System.out.println(CigarOperator.M.consumesReferenceBases()); //True
		System.out.println(CigarOperator.EQ.consumesReferenceBases());
		System.out.println(CigarOperator.X.consumesReferenceBases());
		System.out.println(CigarOperator.N.consumesReferenceBases());
		System.out.println(CigarOperator.D.consumesReferenceBases());
		
		System.out.println(CigarOperator.I.consumesReferenceBases()); //False
		System.out.println(CigarOperator.S.consumesReferenceBases());		
		System.out.println(CigarOperator.P.consumesReferenceBases());
		System.out.println(CigarOperator.H.consumesReferenceBases());
 *
 */
public class AlignedRead {
	
	private byte CHAR_FOR_SOFT_CLIP= (byte)'n';
	private byte[] readbases;
	private byte[] refbases;
	private byte[] methylbases;
	private boolean isForward;
	private boolean isFirst; // First in pair	
	
	/*
	 * S E T T E R S   &   G E T T E R S
	*/

	public byte getCHAR_FOR_SOFT_CLIP() {
		return CHAR_FOR_SOFT_CLIP;
	}
//	public void setCHAR_FOR_SOFT_CLIP(byte char_for_soft_clip) {
//		CHAR_FOR_SOFT_CLIP = char_for_soft_clip;
//	}
	
	public boolean isFirst() {
		return isFirst;
	}
	public byte[] getMethylbases() {
		return methylbases;
	}
	public void setMethylbases(byte[] methylbases) {
		this.methylbases = methylbases;
	}
	public void setFirst(boolean isFirst) {
		this.isFirst = isFirst;
	}
	public byte[] getReadbases() {
		return readbases;
	}
	public void setReadbases(byte[] readSeq) {
		this.readbases = readSeq;
	}

	public byte[] getRefbases() {
		return refbases;
	}
	public void setRefbases(byte[] refSeq) {
		this.refbases = refSeq;
	}

	public boolean isForward() {
		return isForward;
	}
	public void setForward(boolean isForward) {
		this.isForward = isForward;
	}
	
	/* ------------------------
	 *  C O N S T R U C T O R S  
	 --------------------------*/
	
	// IMPORTANT: Set the isForward/isFirst flags 
	// 		      BEFORE running the setRead/Ref/Methy methods!!
	//			  Otherwise you always use false/false (the default boolean)
	
	public AlignedRead(){
		
	};
	
	public AlignedRead(SAMRecord rec){
		this.isForward= !rec.getReadNegativeStrandFlag();
		this.isFirst= setMateFlag(rec);
		this.readbases= this.getAlignedBasesForRead(rec);
	}
	
	public AlignedRead(SAMRecord rec, FastaReference fa){
		this.isForward= !rec.getReadNegativeStrandFlag();
		this.isFirst= setMateFlag(rec);
		this.readbases= this.getAlignedBasesForRead(rec);
		this.refbases= this.getReferenceBasesForRead(rec, fa);
		this.methylbases= this.methylArray(readbases, refbases, isFirst, isForward);
	}

	public AlignedRead(SAMRecord rec, FastaReference fa, byte CHAR_FOR_SOFT_CLIP){
		this.isForward= !rec.getReadNegativeStrandFlag();
		this.isFirst= setMateFlag(rec);
		this.CHAR_FOR_SOFT_CLIP= CHAR_FOR_SOFT_CLIP;
		this.readbases= this.getAlignedBasesForRead(rec);
		this.refbases= this.getReferenceBasesForRead(rec, fa);
		this.methylbases= this.methylArray(readbases, refbases, isFirst, isForward);
	}
	
	/* ------------------------------------------------------------------------
	 *                            M E T H O D S
	 * --------------------------------------------------------------------- */

	private boolean setMateFlag(SAMRecord rec){
		boolean flag;
		if (!rec.getReadPairedFlag()){
			flag= true; // If read is not paired it is treated as first-in-pair
		} else {
			flag= rec.getFirstOfPairFlag();
		}
		return(flag);
	}
	
	private byte[] reverseComplement(byte[] seq){
				
		DNASequence dnaseq= new DNASequence(new String(seq));
		
		return( dnaseq.getReverseComplement().getSequenceAsString().getBytes());

	}
	
	/**
	 * Get the read bases aligned to reference. The output sequence is in the same
	 * orientation as in the bam input. 
	 * @param rec
	 * @return
	 */
	public byte[] getAlignedBasesForRead(SAMRecord rec){
		
		int readPosition= 0;
		byte[] readBases= rec.getReadBases();
		byte[] alignedRead= new byte[readBases.length];
		for (CigarElement cigarElement : rec.getCigar().getCigarElements()){

			CigarOperator cigarOperator= cigarElement.getOperator();
			
			// From sam specs:
			// Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigarOperator.equals(CigarOperator.S)){
				
				for (int i= 0; i < cigarElement.getLength(); i++){
					alignedRead[readPosition]= CHAR_FOR_SOFT_CLIP;
					readPosition++;
				}
				
			} else if (cigarOperator.consumesReadBases()){
				
				for (int i= 0; i < cigarElement.getLength(); i++){
									
					alignedRead[readPosition]= readBases[readPosition];
					readPosition++;
				}
				
			} else {
				continue;
			}
		}
		if (alignedRead.length != readBases.length){
			System.err.println("SEQ length from sam record does not match reconstructed read length.");
			System.exit(1);
		}
				
		return(alignedRead);
	}
	
	/**
	 * Get the reference bases corresponding to each read base. Reference bases appear
	 * in the same orientation as the input bam.
	 * @param rec
	 * @param fa
	 * @return
	 */
	public byte[] getReferenceBasesForRead(SAMRecord rec, FastaReference fa){

		byte[] referenceBases= new byte[rec.getReadLength()];

		List<CigarElement> cigarElements= rec.getCigar().getCigarElements();
		
		int readIdx= 0; 
		
		for (CigarElement el : cigarElements){
			if (el.getOperator().consumesReferenceBases()){

				if (readIdx == 0 && el.getOperator().equals(CigarOperator.D)){
					System.err.println("Warning: Cigar string starts with a deletion!");
				}
				
				for(int i= 0; i < el.getLength(); i++){
										
					int readPos= rec.getReferencePositionAtReadPosition(readIdx+1);
					
					referenceBases[readIdx]= fa.getSubsequenceAt(rec.getReferenceName(), readPos, readPos, true)[0]; 
															
					if (el.getOperator().consumesReadBases()){
						readIdx++;
					}
				}
			} else {
				for(int i= 0; i < el.getLength(); i++){
					if (!el.getOperator().equals(CigarOperator.HARD_CLIP)){
						referenceBases[readIdx]= CHAR_FOR_SOFT_CLIP; 
					}
					if (el.getOperator().consumesReadBases()){
						readIdx++;
					}
				}
			}
				
		}
		for (int i=0; i < referenceBases.length; i++){
			if (referenceBases[i] == 0){
				System.err.println("Reference array contains null byte at index pos: " + i);
				System.err.println(new String(referenceBases));
				System.exit(1);
			}
		}
		
		return(referenceBases);
	}
	
	/**
	 * Look at C on the reference C: 
	 *    - Forward & First in pair
	 *    - Reverse & Second in pair
	 * Look at G on the
	 *    - Reverse & First in pair
	 *    - Forward & Second in pair
	 *    
	 * In other words: Look at C if you have true/true or false/false. Look at G otherwise.
	 * */
	public char nucleotideToMethylation(char base, char ref, boolean isFirst, boolean isForward){
		base= Character.toUpperCase(base);
		ref= Character.toUpperCase(ref);
		if ((isForward && isFirst) || (!isForward && !isFirst)){
			if( ref != 'C') return '.';
			if( ref == 'C' && base == 'C') return 'M';
			if( ref == 'C' && base == 'T') return 'u';
			if( ref == 'C' && (base != 'C' || base != 'T')) return '*';
			System.err.println("Unexpected combination of base and reference");
			System.exit(1);			
		}
		if ((!isForward && isFirst) || (isForward && !isFirst)){
			if( ref != 'G' ) return '.';
			if( ref == 'G' && base == 'G') return 'M';
			if( ref == 'G' && base == 'A') return 'u';
			if( ref == 'G' && (base != 'G' || base != 'A')) return '*';
			System.err.println("Unexpected combination of base and reference");
			System.exit(1);
		}
		System.err.println("Unexpected combination of strand and pairing information");
		System.exit(1);
		return('.');
	}
	
	
	public byte[] methylArray(byte[] read, byte[] ref, boolean isFirst, boolean isForward) {
		
		byte[] methyl= new byte[read.length];
		for (int i=0; i < methyl.length; i++){
			methyl[i]= (byte) nucleotideToMethylation((char)read[i], (char)ref[i], isFirst, isForward);
		}
		return methyl;
	}

	public byte[] methylArray(){
		
		byte[] read= this.getReadbases();
		byte[] ref= this.getRefbases();
		byte[] methyl= methylArray(read, ref, this.isFirst, this.isForward);
		return(methyl);
	}
	
	
	/**
	 * Re-orient the reads, reference, and methylation arrays to have the 
	 * reads appearing as in the original fastq.
	 * I.e. rev comp if mapped to reverse strand. Methyaltion array is simply reversed.
	 */
	public void orientAsInFastq(){
		if (!this.isForward){
			this.readbases= reverseComplement(this.readbases);
			this.refbases= reverseComplement(this.refbases);
			ArrayUtils.reverse(this.methylbases);
		}
	}
	
	/**
	 * Remove from arrays read/ref/methyl the character marked by CHAR_FOR_SOFT_CLIP
	 */
	public void clipAlignment(){
		
		List<Integer> softClippedPositions= this.getSoftClippedPositions(this.getReadbases()); 
		this.readbases= this.softClipArray(softClippedPositions, this.getReadbases());
		this.refbases= this.softClipArray(softClippedPositions, this.getRefbases());
		this.methylbases= this.softClipArray(softClippedPositions, this.getMethylbases());
		
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append(new String(this.getReadbases())); sb.append("\t");
		sb.append(new String(this.getRefbases()));  sb.append("\t");
		sb.append(new String(this.getMethylbases()));  sb.append("\t");
		sb.append(this.isForward); sb.append("\t");
		sb.append(this.isFirst);
		return(sb.toString());
	}
	
	/**
	 * Get positions of soft clipped bases from array of read bases.
	 * @param readBases
	 * @return
	 */
	private List<Integer> getSoftClippedPositions(byte[] readBases){
		
		List<Integer> sofClippedPositions= new ArrayList<Integer>();
		for (int i= 0; i < readBases.length; i++){
			char c= (char)readBases[i];
			if (c == CHAR_FOR_SOFT_CLIP){
				sofClippedPositions.add(i);
			}
		}
		return(sofClippedPositions);
	}
	
	/**
	 * Return a copy of the input array (e.g. readbases) without the positions given in the list 
	 * of integer 
	 * @param softClippedPositions Remove these positions from the array
	 * @param arrayToClip Array to clip
	 * @return
	 */
	private byte[] softClipArray(List<Integer> softClippedPositions, byte[] arrayToClip){
		byte[] clippedArray= new byte[arrayToClip.length - softClippedPositions.size()];
		
		int p= 0;
		for(int i= 0; i < arrayToClip.length; i++){
			if(!softClippedPositions.contains(i)){
				clippedArray[p]= arrayToClip[i];
				p++;
			}
		}
		return(clippedArray);
	}
}