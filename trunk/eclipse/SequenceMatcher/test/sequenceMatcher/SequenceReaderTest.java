package sequenceMatcher;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import org.biojava3.core.sequence.DNASequence;
import org.junit.Test;

import com.google.common.primitives.Chars;

public class SequenceReaderTest {

	@Test
	public void canGetNextSequenceFromFasta() throws FileNotFoundException, IOException{
		
		String fastafile= "test/seqs.fa";
				
		BufferedReader br= new BufferedReader(new FileReader(fastafile));
		String[] faseq= null;
		faseq= SequenceReader.getNextSequence(br);	
		assertEquals(faseq[0], "seq1");
		assertEquals(faseq[1], "ACTGAAATTTCCCCCCTTTCCCTTTACCCTTTCCCTTTA");
		
		faseq= SequenceReader.getNextSequence(br); // seq2
		
		faseq= SequenceReader.getNextSequence(br);
		assertEquals(faseq[0], "seq3");
		assertEquals(faseq[1], "ACTG");
		
		faseq= SequenceReader.getNextSequence(br); // No sequence left to read
		assertNull(faseq);
		br.close();
	}
	
	
	@Test
	public void canOpenFastaFile() throws FileNotFoundException, IOException{

		String fastafile= "test/seqs.fa";
		ArrayList<String[]> fasta= SequenceReader.readFastaToList(fastafile);		
		
		assertEquals(fasta.get(0)[0], "seq1");
		assertEquals(fasta.get(0)[1], "ACTGAAATTTCCCCCCTTTCCCTTTACCCTTTCCCTTTA");
				
	}
	
	@Test
	public void speedDNASequence(){
		String dna= "ACTGNACTGNACTGNACTGNACTGNACTGNACTGN";
				
		// HashMap<Character, Character> iupac= SequenceReader.getIupacAmbiguityDNA(); 

		int i= 0;
		String rc;
		AmbiguityDNACompoundSet ambDNAset= new AmbiguityDNACompoundSet();
		long t0= System.currentTimeMillis();
		while(i < 1000000){
//			rc= SequenceReader.revcomp(dna, iupac);
			rc= new DNASequence(dna, ambDNAset).getReverseComplement().getSequenceAsString();
			i++;
		}
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0);
	}
	
	private void print(Object x){
		System.out.println(x);
	}
}

/*
 	@Test
	public void canReadAllFasta() throws FileNotFoundException, IOException{

		String fastafile= "test/seqs.fa";
		HashMap<String, String[]> farray= SequenceReader.fastaToArray(fastafile);
	
		assertTrue(farray.size() == 2);
	
		int nseq= farray.get("name").length;
		assertTrue(nseq == 3);
		// Name of last sequence is "seq3"
		assertEquals(farray.get("name")[nseq-1], "seq3");
		// Sequence of last entry:
		assertEquals(farray.get("seq")[nseq-1], "ACTG");
		
	}
 */