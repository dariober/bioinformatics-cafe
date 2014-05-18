package markovChain;

import com.google.common.collect.*;

import java.lang.reflect.Array;
import java.util.*;

public class Tester {
		
	public static void main(String[] args){
		
		Model model= new Model(2, new Alphabet().getNucleotides(), 
				"CACACACACACTG", 0.01);
		System.out.println(model);
		
		MarkovSequence mc= new MarkovSequence();
		String mcseq= mc.generateMarkovSequence(model, 1000);
		System.out.println(mcseq);
				
//		Model modelTemplate= new Model();
//		Model model2= new Model(model, "ATATATATA");
//		System.out.println(model2);
//		System.exit(0);
		
		
//		testInitializeModel();
//		timeitGenerateModel();
//		timeitInitializeChain();
//		timeitGetNextChar();
//		timeitGenerateMarkovSequence();
	}
}
