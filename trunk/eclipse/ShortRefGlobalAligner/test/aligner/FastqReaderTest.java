package aligner;

import static org.junit.Assert.*;

import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class FastqReaderTest {
	
	@Test
	public void canOpenFastqFile() throws FileNotFoundException, IOException{
		BufferedReader br= FastqReader.openFastq("test_data/fk018.fq.gz");
		String header= br.readLine();
		assertTrue(header.startsWith("@"));
		br.close();
		
		br= FastqReader.openFastq("test_data/fk018.fq");
		header= br.readLine();
		assertTrue(header.startsWith("@"));
		br.close();
	}
	
	@Test
	public void readFirstFourLines() throws FileNotFoundException, IOException{
		
		/*
		@M00886:24:000000000-A9PNT:1:1101:19356:1061 1:N:0:3
		NAGTTTCGATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATTAGAGGACCCGTGTCAGATATATACATCCGAT
  		+
        #8ACCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGG
		*/
		
		BufferedReader br= FastqReader.openFastq("test_data/fk018.fq.gz");
		String[] fqread= FastqReader.getNextRead(br);
		
		assertEquals("M00886:24:000000000-A9PNT:1:1101:19356:1061 1:N:0:3", fqread[0]);
		assertEquals("NAGTTTCGATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATTAGAGGACCCGTGTCAGATATATACATCCGAT", fqread[1]);
		assertEquals("", fqread[2]);
		assertEquals("#8ACCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGG", fqread[3]);
		
		fqread= FastqReader.getNextRead(br);
		assertEquals("M00886:24:000000000-A9PNT:1:1101:13784:1062 1:N:0:3", fqread[0]);
		br.close();
	}
	
	@Test
	public void readToEndOfFile() throws FileNotFoundException, IOException{
		
		BufferedReader br= FastqReader.openFastq("test_data/fk018.fq.gz");
		String[] fqread= FastqReader.getNextRead(br);
		
		List<String> reads= new ArrayList<String>();
		while(fqread != null){
			reads.add(fqread[0]);
			fqread= FastqReader.getNextRead(br);
		}
		assertEquals(250, reads.size()); // Number of reads in fastq
		assertEquals("M00886:24:000000000-A9PNT:1:1101:12195:1559 1:N:0:3", reads.get(reads.size()-1)); // Name of last read 
		br.close();
	}
}
