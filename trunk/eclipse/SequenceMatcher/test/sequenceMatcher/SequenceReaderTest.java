package sequenceMatcher;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.DNASequence;
import org.junit.Test;


public class SequenceReaderTest {

	@Test
	public void canOpenFastaFile() throws FileNotFoundException, IOException{

		String fastafile= "test/seqs.fa";
		LinkedHashMap<String, DNASequence> fasta= SequenceReader.ReadFasta(fastafile);			
		assertEquals(fasta.get("seq3").getSequenceAsString(), "ACTG");
		
	}

	@Test
	public void canConvertFastaHashMapToHashMapOfStringArrays() throws FileNotFoundException, IOException{

		String fastafile= "test/seqs.fa";
		LinkedHashMap<String, DNASequence> fasta= SequenceReader.ReadFasta(fastafile);
		HashMap<String, String[]> fastamap= SequenceReader.DNASequenceMapToArrays(fasta);

		// All elements read:
		assertTrue(fastamap.get("name").length == fasta.size());
		assertTrue(fastamap.get("seq").length == fasta.size());

		// Get right elements:
		assertEquals(fastamap.get("name")[0], "seq1");
		assertEquals(fastamap.get("seq")[2], "ACTG");
	}
}
