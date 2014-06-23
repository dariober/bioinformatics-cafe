package bisReadBias;

import java.util.List;

import org.biojava3.core.sequence.DNASequence;

// import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * @author berald01
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
	
	private byte CHAR_FOR_SOFT_CLIP= (byte)'N';
	private byte[] readbases;
	private byte[] refbases;
	private byte[] methylbases;
	private boolean isForward;
	private boolean isFirst; // First in pair	
	
	/*
	 * S E T T E R S   &   G E T T E R S
	*/
	
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
		
	public AlignedRead(){
		
	};
	
	public AlignedRead(SAMRecord rec){
		this.readbases= this.getAlignedBasesForRead(rec);
		this.isForward= !rec.getReadNegativeStrandFlag();
		this.isFirst= setMateFlag(rec);
	}
	
	public AlignedRead(SAMRecord rec, FastaReference fa){
		this.readbases= this.getAlignedBasesForRead(rec);
		this.refbases= this.getReferenceBasesForRead(rec, fa);
		this.isFirst= setMateFlag(rec);
		this.methylbases= this.methylArray(readbases, refbases);
	}
	
	/* ------------------------------------------------------------------------
	 *                            M E T H O D S
	 * --------------------------------------------------------------------- */

	private boolean setMateFlag(SAMRecord rec){
		boolean flag;
		if (!rec.getReadPairedFlag()){
			flag= true;
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
	 * order as in the fastq read. 
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
				
		if (rec.getReadNegativeStrandFlag()){
			alignedRead= reverseComplement(alignedRead);
		}
		
		return(alignedRead);
	}
	
	/**
	 * Get the reference bases corresponding to each read base. Reference bases appear
	 * in the same order as the fastq read.
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
					referenceBases[readIdx]= CHAR_FOR_SOFT_CLIP; 
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
		
		if (rec.getReadNegativeStrandFlag()){
			referenceBases= reverseComplement(referenceBases);
		}
		
		return(referenceBases);
	}
	
	public byte[] methylArray(byte[] read, byte[] ref) {
		
		byte[] methyl= new byte[read.length];
		for (int i=0; i < methyl.length; i++){

			char xread= (char) Character.toUpperCase(read[i]);
			char xref= (char) Character.toUpperCase(ref[i]);
			char xmeth;
			Character.toUpperCase(xref);
			if (xref == 'C'){
				if (xread == 'C'){
					xmeth= 'M';
				} else if (xread == 'T'){
					xmeth= 'u';
				} else {
					// Mismatch
					xmeth= '*';
				}
			} else {
				// Any base other than C matching or not
				xmeth= '.';
			}
			methyl[i]= (byte)xmeth;
		}
		return(methyl);
		
	}
	public byte[] methylArray(){
		
		byte[] read= this.getReadbases();
		byte[] ref= this.getRefbases();
		byte[] methyl= methylArray(read, ref);
		return(methyl);
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
}