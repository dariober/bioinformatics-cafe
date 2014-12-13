package sequenceMatcher;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashMap;

import org.junit.Test;


public class DistanceTest {

	@Test
	public void hammingOnEqualSeqs(){
	
		String seqA= "ACTG";
		String seqB= "ACTG";
		int d= Distance.hammingDist(seqA, seqB, 0);
		assertTrue(d == 0);
		
		d= Distance.hammingDist(seqA, seqB, 1);
		assertTrue(d == 0);

		d= Distance.hammingDist(seqA, seqB, -1);
		assertTrue(d == 0);
		
	}
	
	@Test
	public void hammingOnDifferentSeqs(){
		
		String seqA= "ACTGNN";
		String seqB= "ACTG";
		int d= Distance.hammingDist(seqA, seqB, -1);
		assertTrue(d == 2);

		d= Distance.hammingDist(seqA, seqB, 0);
		assertTrue(d == -1);

		d= Distance.hammingDist(seqA, seqB, 1);
		assertTrue(d == -1);
		
		d= Distance.hammingDist(seqA, seqB, 2);
		assertTrue(d == 2);

		d= Distance.hammingDist(seqA, seqB, 3);
		assertTrue(d == 2);
		
	}

	@Test
	public void canGetNMtagFromAln(){
		String alnRead;
		String alnRef;
		int nm= -1;
		
		alnRead= "ACTG";
		alnRef=  "ACTG";
		nm= Distance.getNMtagFromAln(alnRead, alnRef);
		assertEquals(0, nm);

		alnRead= "ACTGA";
		alnRef=  "ACTGN";
		nm= Distance.getNMtagFromAln(alnRead, alnRef);
		assertEquals(1, nm);

		alnRead= "AACTG-A";
		alnRef=  "A-CTGTN";
		nm= Distance.getNMtagFromAln(alnRead, alnRef);
		assertEquals(3, nm);

		alnRead= "AACTG--";
		alnRef=  "AACTGNN";
		nm= Distance.getNMtagFromAln(alnRead, alnRef);
		assertEquals(0, nm);
		
		alnRead= "---AACTG--";
		alnRef=  "NNNACCTGNN";
		nm= Distance.getNMtagFromAln(alnRead, alnRef);
		assertEquals(1, nm);
	
	}
	
}
