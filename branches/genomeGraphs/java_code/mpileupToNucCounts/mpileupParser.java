
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.String;
import java.util.*;

/*
*  _MEMO_
* Compile java
* ------------
* cd /to/mpileupParser.java
* javac mpileupParser.java
*
* 
* Create jar executable with
* --------------------------
* cd /dir/to/with-all-*.class-files
* jar cmf manifest.txt mpileupToNucCounts.jar *.class
*
* Where manifest.txt contains one line with (must end with newline):
* ------
* Main-Class: mpileupParser
*
* ------
* mpileupParser is the class with function `public static void main`
*/


class Pile {

	public String extractLineCoords(String[] pileupLine){
		
		/*
		*   Connvert coordinates and refbase to String formatted as (partial) 
		*   python dict
		*/
		
		String pydict= "{'chrom': '" + pileupLine[0] + "', 'pos': " + pileupLine[1] + ", 'base': '" + pileupLine[2] + "', ";
		return(pydict);
	}

	public Map<String, Integer> condenseBasesToHashMapCount(String bases, String[] nucs, String refbase){
		/*
		*   Condense one string of bases from mpileup (not the whole line) to base
		*   counts.
		*   
		*   bases:
		*   	String of bases from mpileup line
		*   nucs:
		*   	List of nucleotides to call
		*   refbase:
		*   	Reference base extracted from 3rd column of mpileup
		*
		*	Return:
		*		HashMap of base calls:
		*		{T=0, g=0, t=0, G=0, c=0, A=0, C=0, a=0, n=0, N=0, .=0, ,=0}
		*
		*   MEMO: mpileip line looks like this:
		*  
		*   chr7	3282070	N	0	*	*	0	*	*	2	^~G^~G	AA
		*	chr7	3282071	N	0	*	*	0	*	*	2	GG	BA
		*
		*/

		Map<String, Integer>  callDict = new HashMap<String, Integer>();

		for(String k : nucs){
			callDict.put(k, 0);
		}

		boolean skip= false;

		for (int i = 0, n = bases.length(); i < n; i++) {
			String b = Character.toString(bases.charAt(i));
			if(b  == "^"){
				skip= true;
			} else if(skip){
				skip= false;
			} else if(callDict.containsKey(b)){
				callDict.put(b, callDict.get(b) + 1);
			}
			else {
				continue;
			}
		}

		callDict.put(refbase,               callDict.get(refbase) +               callDict.get("."));
		callDict.put(refbase.toLowerCase(), callDict.get(refbase.toLowerCase()) + callDict.get(","));

		return(callDict);
	}

    public List<Map> extractNucCountsFromPileupLine (String[] pileupLine) {

		/*
		*   Get a whole line of input from samtools mpileup and return a list of
		*   hashes, one hash for each bam file
		*/
	
		String[] nucs= {"A", "a", "C", "c", "G", "g", "T", "t", "N", "n", ".", ","};
		
		String refbase= pileupLine[2];
		
		String bases= "";
		Pile d = new Pile();
		List<Map> listOfNucCounts = new ArrayList<Map>(); // List with one hashmap per bamfile
		for(int i = 4; i < pileupLine.length; i = i+3) {
			bases= pileupLine[i];
			Map <String, Integer> dict = d.condenseBasesToHashMapCount(bases, nucs, refbase);
			listOfNucCounts.add(dict);
		}
		return(listOfNucCounts);
    }
	
	public String countMapToPydict(List<Map> listOfNucCounts){
		
		/*
		*  Convert list of MapHashes containing nuc counts to String formatted
		*  as (partial) python dict
		*/
		
		String[] forwardNucs= {"A", "C", "G", "T", "N"};
		String[] revNucs=     {"a", "c", "g", "t", "n"};
		String[] outnucs=     {"A", "a", "C", "c", "G", "g", "T", "t", "N", "n"}; // MEMO: Add "Z" and "z" 

		String pydict = "";
		int i= 0;
		for(Map nucCountMap : listOfNucCounts){
			int Z= 0;
			for(String k : forwardNucs){
				int x = (Integer) nucCountMap.get(k);
				Z += x;
			}
			int z= 0;
			for(String k : revNucs){
				int x = (Integer) nucCountMap.get(k);
				z += x;
			}
			String innerDict= i + ": {";
			for(String x : outnucs){
				innerDict += "'" + x + "': " + nucCountMap.get(x) + ", "; // Memo: No need to strip last comma {1:0, }
			}
			innerDict += "'Z': " + Z + ", 'z': " + z + "}, ";
			i += 1;
			pydict += innerDict;
		}
		pydict += "}";
		return(pydict);
	}
}

public class mpileupParser{
 
	public static void main (String args[]) {
 
		try{
			BufferedReader br = 
						  new BufferedReader(new InputStreamReader(System.in));
			
			String input;	
			Pile d = new Pile();			
			while((input=br.readLine())!=null){
				String[] pileupLine= input.split("\t");
				String coords= d.extractLineCoords(pileupLine);
				List<Map> listOfNucCounts = d.extractNucCountsFromPileupLine(pileupLine);
				String counts = d.countMapToPydict(listOfNucCounts);
				System.out.println(coords + counts);
			}		
		}
		catch(IOException io) {
			io.printStackTrace();
		}
	}
}
