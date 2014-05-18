package regexEnrichment;

import java.util.*;

/**
 * Class to hold a contig sequence divided into chunks of fixed size.
 * Store also the start and end of each chunk.
 * @author berald01
 *
 */
public class ContigWindows {

	public ContigWindows() {
		this.sequenceName = null;
		this.sequence = null;
	}
	
	public ContigWindows(String sequenceName, String sequence) {
		this.sequenceName = sequenceName;
		this.sequence = sequence;
	}
	
	private String sequenceName;
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	private String sequence;	
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	private int precision;
	public String getPrecision() {
		return sequence;
	}
	public void setPrecision(int precision) {
		this.precision = precision;
	}
	
	private List<Window> windowCoords= new ArrayList<Window>();
	
	public List<Window> setWindowCoords(int windowSize, int step){

		int firstEnd= (windowSize > sequence.length()) ? sequence.length() : windowSize;

		Window window= new Window();		
		window.setChrom(sequenceName);
		window.setStart(0);
		window.setEnd(firstEnd);
		window.setSequence(sequence.substring(0, firstEnd));
		window.setNucleotideFreq(precision);
		windowCoords.add(window);
		int start= step; // windowSize;
		while (start < sequence.length()){
			int end= ((start + windowSize) > sequence.length()) ? sequence.length() : (start + windowSize);
			window= new Window();
			window.setChrom(sequenceName);
			window.setStart(start);
			window.setEnd(end);
			window.setSequence(sequence.substring(start, end));
			window.setNucleotideFreq(precision);
			windowCoords.add(window);
			start += step; // windowSize;
		}
		return (windowCoords);
	}
	
	public List<Window> getWindowCoords(){
		return (windowCoords);
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();;
		for(Window w : this.getWindowCoords()){
			sb.append(w.toString());
			sb.append("\n");
		}
		return(sb.toString().trim());
	}
}


