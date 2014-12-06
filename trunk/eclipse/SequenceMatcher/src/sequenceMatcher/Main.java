package sequenceMatcher;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Set;

import org.biojava3.core.sequence.DNASequence;

public class Main {

	public static void main(String[] args) throws IOException {
		
		int NMAX= 100;
		
		LinkedHashMap<String, DNASequence> set1= SequenceReader.ReadFasta("test/seqs.fa");
		LinkedHashMap<String, DNASequence> set2= SequenceReader.ReadFasta("test/seqs.fa");
		
		HashMap<String, String[]> array1= SequenceReader.DNASequenceMapToArrays(set1);
		HashMap<String, String[]> array2= SequenceReader.DNASequenceMapToArrays(set1);
		
		String s1= null;
		String s2= null;
		for(int i=0; i < set1.size(); i++){
			s1= array1.get("seq")[i];
			for(int j=0; j < set2.size(); j++){
				s2= array2.get("seq")[j];
				int d= Distance.levenDist(s1, s2, NMAX);
				if(d <= NMAX){
					System.out.println(s1 + " " + s2 + " " + d);
				}
			}
		}

	}

}
