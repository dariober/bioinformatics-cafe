package sequenceMatcher;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.zip.GZIPInputStream;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.DNASequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;

public class SequenceReader {

	/**
	 * Read entire fasta file and put it in LinkedHashMap.
	 * @param fastafile
	 * @return
	 * @throws IOException
	 */
	private static LinkedHashMap<String, DNASequence> readFasta(String fastafile) throws IOException{
		
		InputStream inStream= Opener.open(fastafile);
	
		GenericFastaHeaderParser<DNASequence, NucleotideCompound> genericFastaHeaderParser= 
				new GenericFastaHeaderParser<DNASequence, NucleotideCompound>();
		DNASequenceCreator dnaSequenceCreator= 
				new DNASequenceCreator(AmbiguityDNACompoundSet.getDNACompoundSet());
	
		FastaReader<DNASequence, NucleotideCompound> fastaReader=
				new FastaReader<DNASequence, NucleotideCompound>(
						inStream, genericFastaHeaderParser, dnaSequenceCreator); 
		LinkedHashMap<String, DNASequence> fastaseq = fastaReader.process();
			
		return(fastaseq);
	}

	/**
	 * Convert a LinkedHashMap containing fasta sequences to a List of String arrays.
	 * Each array has length 3 and contains a fasta sequence: {name, sequence, sequence revcomp}.
	 * The LinkedHashMap is typically obtained by calling ReadFasta(fastafile).
	 * @param fastaMap
	 * @return
	 */
	private static ArrayList<String[]> fastaLinkedHashMapToList(
			LinkedHashMap<String, DNASequence> fastaMap){
		
		ArrayList<String[]> fastaList= new ArrayList<String[]>();
		AmbiguityDNACompoundSet ambDNAset= new AmbiguityDNACompoundSet();
		
		for(String seqname : fastaMap.keySet()){						
			String[] faseq= new String[3];
			faseq[0]= seqname;
			faseq[1]= fastaMap.get(seqname).getSequenceAsString().toUpperCase();
			faseq[2]= new DNASequence(faseq[1], ambDNAset).getReverseComplement().getSequenceAsString();
			fastaList.add(faseq);
		}
		
		return fastaList;
	}
	
	/**
	 * Read a fasta file and put each sequence (name & sequence pair) in a list of
	 * String arrays. Each array has length 2: {name, sequence}.
	 * @param fastafile
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<String[]> readFastaToList(String fastafile) throws IOException{
		ArrayList<String[]> fastaList= fastaLinkedHashMapToList(readFasta(fastafile));
		return fastaList;
	}
	
	public static BufferedReader openFastq(String fastq) throws FileNotFoundException, IOException{
		
		BufferedReader br;
		if(fastq.endsWith(".gz")){
		    br=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
		} else {
		   br=new BufferedReader(new FileReader(fastq));
		}
		return(br);
	}

	/**
	 * @param br
	 * @return String array of length 4 with 
	 * 1) read name (w/o leading @)
	 * 2) Sequence
	 * 3) Comment (w/o leading +)
	 * 4) Quality
	 * @throws IOException
	 */
	public static String[] getNextRead(BufferedReader br) throws IOException {
		String[] fqread= new String[4];
		
		// Name
		String name= br.readLine();
		if (name == null){
			return(null);
		}
		fqread[0]= name.substring(1, name.length());
		// Sequence
		fqread[1]= br.readLine();
		// Comment
		String comment= br.readLine();
		fqread[2]= comment.substring(1, comment.length());
		// Quality
		fqread[3]= br.readLine();	
		return(fqread);
	}

	/**
	 * Read next sequence from FASTA file and put it in a String array 
	 * of length two:
	 * String[0]: Name
	 * String[1]: Sequence
	 * @param br BufferedReader connected to the fasta file to read.
	 * @return
	 * @throws IOException
	 */
	public static String[] getNextSequence(BufferedReader br) throws IOException {
		
		String[] fastaseq= new String[2];
		
		int BUFFER_SIZE = 8192;	
		String line= br.readLine();

		// Name
		if(line == null){
			return null;
		} else if (line.startsWith(">")){
			fastaseq[0]= line.replaceFirst(">", "");
		} else {
			System.err.println(line);
			System.err.println("Invalid sequence name or format");
			System.exit(1);
		}
		StringBuilder sb= new StringBuilder();
		while(true){
			br.mark(BUFFER_SIZE);
			line= br.readLine();
			if(line == null || line.startsWith(">")){
				break;
			} else {
				sb.append(line);
			}
		}		
		String sequence= sb.toString().toUpperCase();
		fastaseq[1]= sequence;
		br.reset();
		return fastaseq;		
	}
	
	public static String revcomp(String dna, HashMap<Character, Character> iupac){
		char[] f= dna.toCharArray();
		char[] rc= new char[f.length];

		for(int i= 0; i < f.length; i++){
			rc[i]= iupac.get(f[i]);
		}
		return(Arrays.toString(rc));
	}
	
}

/*	
private HashMap<Character, Character> iupacAmbiguityDNA= new HashMap<Character, Character>();
	
public void setIupacAmbiguityDNA(HashMap<Character, Character> iupacAmbiguityDNA) {
	this.iupacAmbiguityDNA = iupacAmbiguityDNA;
}

public static HashMap<Character, Character> getIupacAmbiguityDNA() {

HashMap<Character, Character> iupacAmbiguityDNA= new HashMap<Character, Character>();
iupacAmbiguityDNA.put('A', 'T');
iupacAmbiguityDNA.put('G', 'C');
iupacAmbiguityDNA.put('C', 'G');
iupacAmbiguityDNA.put('T', 'A');
iupacAmbiguityDNA.put('Y', 'R');
iupacAmbiguityDNA.put('R', 'Y');
iupacAmbiguityDNA.put('W', 'W');
iupacAmbiguityDNA.put('S', 'S');
iupacAmbiguityDNA.put('K', 'M');
iupacAmbiguityDNA.put('M', 'K');
iupacAmbiguityDNA.put('D', 'H');
iupacAmbiguityDNA.put('V', 'B');
iupacAmbiguityDNA.put('H', 'D');
iupacAmbiguityDNA.put('B', 'V');
iupacAmbiguityDNA.put('X', 'X');
iupacAmbiguityDNA.put('N', 'N');

return iupacAmbiguityDNA;
}
*/
