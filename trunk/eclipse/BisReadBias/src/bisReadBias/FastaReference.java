package bisReadBias;

import java.io.File;
import java.util.*;
import java.util.List;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;

public class FastaReference {

	private HashMap<String, byte[]> reference= new HashMap<String, byte[]>();
	
	/* SETTERS AND GETTERS 
	 * -----------------------------------------------------------------------*/
	
	public HashMap<String, byte[]> getReference() {
		return reference;
	}

	public void setReference(HashMap<String, byte[]> reference) {
		this.reference = reference;
	}

	/* CONSTRUCTORS
	 * ----------------------------------------------------------------------*/
	
	public FastaReference(File fastafile) {
		
		try {
			HashMap<String, DNASequence> referenceDNAseq = 
					FastaReaderHelper.readFastaDNASequence(fastafile);
			for (String chrom : referenceDNAseq.keySet()){
				this.reference.put(chrom, DNASequenceToByteArray(referenceDNAseq.get(chrom)));
			}
		
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/* METHODS
	 * ----------------------------------------------------------------------*/
	
	private byte[] DNASequenceToByteArray(DNASequence dnaseq){
		
		List<NucleotideCompound> nucCompList= dnaseq.getAsList();
		
		byte[] seqarray= new byte[nucCompList.size()];
		
		for (int i=0; i < nucCompList.size(); i++){
			NucleotideCompound nc= nucCompList.get(i);
			seqarray[i]= nc.getBase().getBytes()[0];
		}
		return(seqarray);
	}
	
	
	/**
	 * Extract subsequence from ref. Similar to picard IndexedFastaFile.getSubsequenceAt()
	 * NB: "from" and "to" params are 1-based. 
	 *  * The first base is accessed with from=1, to=1
	 *  * The first two bases with from=1, to=2 etc. 
	 *  * If ref length= 100, the last base is accessed with from=100, to= 100; 
	 * @param refname
	 * @param from
	 * @param to
	 * @param allowRightPadding If argument passed to "to" goes beyond the end of the sequence, 
	 * allow right padding with 'N';
	 * @return
	 */
	public byte[] getSubsequenceAt(String refname, int from, int to, boolean allowRightPadding) {
		
		if(from < 1 || from > to){
			System.err.println("Invalid from to args");
			System.exit(1);
		}

		byte[] subsequence= new byte[(to - from) + 1];
		
		byte[] refsequence= this.reference.get(refname);

		int subidx= 0;
		for(int i= subidx; i < subsequence.length; i++){

			int positionOnReference= i+from-1;
			
			byte b;
			
			if(allowRightPadding && positionOnReference >= refsequence.length){
				b= (byte)'N';
			} else {
				b= refsequence[positionOnReference];
			}
			
			subsequence[subidx]= b;
			subidx++;
		}
		return(subsequence);
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		for(String chrom : this.reference.keySet()){
			
			sb.append(chrom + ": ");
			String seq= new String(this.reference.get(chrom));
			sb.append(seq.length() + "bp: ");
			int trimto= (seq.length() > 50) ? 50 : seq.length();
			sb.append(seq.substring(0, trimto) + "\n");
			
		}
		return(sb.toString());
	}

}
