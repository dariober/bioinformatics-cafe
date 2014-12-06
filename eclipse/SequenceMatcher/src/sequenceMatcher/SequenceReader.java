package sequenceMatcher;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
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

	public static LinkedHashMap<String, DNASequence> ReadFasta(String fastafile) throws IOException{
	
		FileInputStream inStream = new FileInputStream( fastafile );
	
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
	 * Convert LinkedHashMap to a pair of String arrays contained in a HashMap.
	 * Format is: 
	 * {"name": ['seq1', 'seq2', ...], "seq": ["ACTG", "AAA", ...]}
	 * @param fastaseq
	 * @return HashMap with key/value pairs: 
	 *     "name": Array of sequence names
	 *     "seq": Array of sequences as strings.
	 */
	public static HashMap<String, String[]> DNASequenceMapToArrays(LinkedHashMap<String, DNASequence> fastaseq){
		
		// Array of sequence names:
		String[] names= fastaseq.keySet().toArray(new String[fastaseq.size()]);

		String[] seqs= new String[fastaseq.size()];
		int i= 0;
		for(String k : fastaseq.keySet()){
			seqs[i]= fastaseq.get(k).getSequenceAsString();
			i++;
		}
		HashMap<String, String[]> fastaMap= new HashMap<String, String[]>();
		fastaMap.put("name", names);
		fastaMap.put("seq", seqs);
		
		return fastaMap;
		
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
	
}
