package bisReadBias;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Joiner;
import com.google.common.collect.*;

/**
 	Data struct
	For read 1 and read2: List of HashMultisets. 
	Each item in List is a position along the read. 
	Each HashMultiset has counts for Methylate, Unmethylated, non-C, mismatch at C.
	E.g.
	Read1 - a List where each element is a position:
	pos 0 = [M:10, u:5, *:1]
	pos 1 = [M:10, u:6, *:1, .:5]
	 	.
	 	N
 * 
 * @author berald01
 *
 */
public class ReadProfile {
	
	private char[] METHYL_CODES= new char[] {'M', 'u', '*', '.'};
	
	private List<HashMultiset<Character>> read1= new ArrayList<HashMultiset<Character>>();
	private List<HashMultiset<Character>> read2= new ArrayList<HashMultiset<Character>>();

	/*
	 * SETTERS AND GETTERS
	 */
	public List<HashMultiset<Character>> getRead1() {
		return read1;
	}

	public void setRead1(List<HashMultiset<Character>> read1) {
		this.read1 = read1;
	}

	public List<HashMultiset<Character>> getRead2() {
		return read2;
	}

	public void setRead2(List<HashMultiset<Character>> read2) {
		this.read2 = read2;
	}

	/**
	 * Increment each position in read 1 or read 2 with the counts in
	 * methylbases
	 * @param Characters
	 * @param isRead2 
	 */
	public void addReadMethyl(byte[] methylbases, boolean isRead1) {
		for (int i=0; i < methylbases.length; i++){
			Character xmeth= (char)methylbases[i];
			if (isRead1){
				if (read1.size() <= i){
					HashMultiset<Character> position= HashMultiset.create();
					read1.add(position);
				}
				read1.get(i).add(xmeth);
			} else {
				if (read2.size() <= i){
					HashMultiset<Character> position= HashMultiset.create();
					read2.add(position);
				}
				read2.get(i).add(xmeth);
			}
		}
	}

	public void addReadMethyl(AlignedRead aln) {
		addReadMethyl(aln.getMethylbases(), aln.isFirst());
	}	
	
	private String printer(HashMultiset<Character> position){
		List<Integer> line= new ArrayList<Integer>();
		for (char x : METHYL_CODES){
			line.add(position.count(x));
		}
		String outline= Joiner.on("\t").join(line);
		return(outline);
	}
	
	public String toString(){

		List<String> outList= new ArrayList<String>();

		// Header order should match METHYL_CODES
		String H= "read\tposition\tcnt_met\tcnt_unmet\tcnt_mismatch\tcnt_noncytosine";
		
		outList.add(H);

		// This complicated list is to iterarate thorugh read1 and read2
		List<ArrayList<HashMultiset<Character>>> reads= new 
					ArrayList<ArrayList<HashMultiset<Character>>>();
		reads.add((ArrayList<HashMultiset<Character>>) this.read1);
		reads.add((ArrayList<HashMultiset<Character>>) this.read2);
		
		int nread= 1;
		for (ArrayList<HashMultiset<Character>> read : reads){
			int pos= 1;
			for(HashMultiset<Character> position : read){
				outList.add(nread + "\t" + pos + "\t" + printer(position));
				pos++;
			}
			nread++;
		}
		
		String outHeader= Joiner.on("\n").join(outList);

		return(outHeader);
	}
}
