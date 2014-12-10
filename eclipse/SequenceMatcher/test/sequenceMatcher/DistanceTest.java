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

}
