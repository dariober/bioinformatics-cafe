package regexEnrichment;

import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.*;
import java.io.*;


// For each nucleotide composition (Map<String, Double>) 
// store the number of hits for each shuffled sequence
public class ProbsLookUpTable implements Serializable {

	private int nrand;
	private String pattern;
	private int window_size;
	private int precision;
	public Map<Map<String, Float>, TreeMultiset<Integer>> 
		probsLookUpTable = new HashMap<Map<String, Float>, TreeMultiset<Integer>>();

	// Keep the actual max num of obs hits. for each combination of obs frequencies.
	// If you encounter
	// a window with more hits than this many hits, you need to run
	// more simulations (unless nrand is already saturated).
	public Map<Map<String, Float>, Integer> 
		maxHitsLookUpTable = new HashMap<Map<String, Float>, Integer>();
	
	// See http://stackoverflow.com/questions/2288937/
	// what-does-it-mean-the-serializable-class-does-not-declare-a-static-final-serial
	private static final long serialVersionUID = 1L; 
	
	ProbsLookUpTable(){}
	
	ProbsLookUpTable (int nrand, String pattern, int window_size, int precision, int maxObsHits){
		this.nrand= nrand;
		this.pattern= pattern;
		this.window_size=window_size;
		this.precision=precision;
	}
	
	public int getNrand() {
		return nrand;
	}
	public void setNrand(int nrand) {
		this.nrand = nrand;
	}
	public String getPattern() {
		return pattern;
	}
	public void setPattern(String pattern) {
		this.pattern = pattern;
	}
	public int getWindow_size() {
		return window_size;
	}
	public void setWindow_size(int window_size) {
		this.window_size = window_size;
	}
	public int getPrecision() {
		return precision;
	}
	public void setPrecision(int precision) {
		this.precision = precision;
	}

	/**
	 * Return true if the number of observed hits for the given nucleotideFreq
	 * is greater than the number of hits encountered so for.
	 * If true, you might need to run more randomisations.
	 * @param obsHits
	 * @param nucleotideFreq
	 * @return
	 */
	public boolean obsHitsGreaterThanMax(int obsHits, Map<String, Float> nucleotideFreq){
		if (!this.maxHitsLookUpTable.containsKey(nucleotideFreq)){
			return(false);
		}
		else {
			int currentMaxHits= this.maxHitsLookUpTable.get(nucleotideFreq);
			return(obsHits > currentMaxHits);
		}
	}
	
	public void writeMap(String filename){
		FileOutputStream fos = null;
		ObjectOutputStream out = null;
		try {
			fos = new FileOutputStream(filename);
			out = new ObjectOutputStream(fos);
			out.writeObject(this);
			out.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	public ProbsLookUpTable readMap(String filename){
		FileInputStream fis = null;
		ObjectInputStream in = null;
		ProbsLookUpTable p= new ProbsLookUpTable();
		try {
			fis = new FileInputStream(filename);
			in = new ObjectInputStream(fis);
			p = (ProbsLookUpTable) in.readObject();
			in.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return(p);
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		
		sb.append("nrand: "); sb.append(nrand);
		sb.append("; pattern: "); sb.append(pattern);
		sb.append("; window_size: "); sb.append(window_size);
		sb.append("; precision: "); sb.append(precision);
		sb.append("; probsLookUpTable: "); sb.append(probsLookUpTable);
		sb.append("; maxHitsLookUpTable: "); sb.append(maxHitsLookUpTable);
		return(sb.toString());
	}
}