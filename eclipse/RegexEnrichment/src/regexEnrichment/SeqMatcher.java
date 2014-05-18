package regexEnrichment;

import java.util.regex.*;

/**
 * Each instance of this class contains a match of the pattern regex to the 
 * target sequence. 
 * 
 * @author berald01
 *
 */
public class SeqMatcher {

	private String group; // Matched group
	private boolean strand; // Strand 
	private int start; // Start/end of the match 
	private int end; 
	
	SeqMatcher(){
		this.strand= true;
	}

	SeqMatcher(String sequence, Matcher m, boolean mstrand){
		this.group= m.group();
		this.strand= mstrand;
		if (this.strand){
			this.start= m.start();
			this.end= m.end();
			this.group= m.group();
		} else {
			this.start= sequence.length() - (m.start() + m.group().length());
			this.end= sequence.length() - m.start();
			this.group= Tools.reverseComplement(m.group());
		}
	}
	
	public String getGroup() {
		return group;
	}

	public boolean isStrand() {
		return strand;
	}

	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
		
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append(this.start); sb.append("\t");
		sb.append(this.end); sb.append("\t");
		sb.append( (this.strand) ? "+" : "-" ); sb.append("\t");
		sb.append(this.group);
		return(sb.toString());
	}
}
