package regexEnrichment;

import java.util.HashMap;

/**
 * Hold data relative to a genomic window, similar to a bed interval 
 * @author berald01
 *
 */
public class Window {

	private String chrom;
	private int start;
	private int end;
	private String sequence;
	private HashMap<String, Float> nucleotideFreq= new HashMap<String, Float>(); 
	
	public String getChrom() {
		return chrom;
	}
	public void setChrom(String contig) {
		this.chrom = contig;
	}

	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		if (start < 0){
			System.out.println("Invalid start: < 0");
			System.exit(1);
		}
		this.start = start;
	}

	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		if (end < 0 || end < start){
			System.out.println("Invalid end: < 0 or >start");
			System.exit(1);
		}
		this.end = end;
	}
	
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		if (sequence.length() != (end - start)){
			System.out.println("Sequence length does not match `end - start`");
			System.exit(1);
		}
		this.sequence = sequence;
	}

	public HashMap<String, Float> getNucleotideFreq() {
		return nucleotideFreq;
	}
	public void setNucleotideFreq(int precision) {
		this.nucleotideFreq= Tools.nucleotideFrequencies(sequence, precision);
	}	
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append("Contig: "); sb.append(chrom);
		sb.append("; Start: "); sb.append(start);
		sb.append("; End: "); sb.append(end);
		sb.append("; Nt. freq.: "); sb.append(nucleotideFreq);
		sb.append("; Sequence: "); sb.append( 
				(sequence != null && sequence.length() > 10 ) ? 
						sequence.substring(0, 10) + "..." : 
							sequence);		
		return(sb.toString());
	}
}
