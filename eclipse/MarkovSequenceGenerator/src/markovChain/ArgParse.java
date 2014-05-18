package markovChain;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import java.util.*;

public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("MarkovChainGenerator")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
						+ "Generate random sequences from input sequence using"
						+ "Markov chain model"
						+ "\n\n"
						+ "EXAMPLE:\n"

						+ "\n\n"
						+ "OUTPUT:\n"

						+ "");
		parser.addArgument("-f", "--fasta")
			.type(String.class)
			.required(true)
			.help("Input fasta file to generate sequences from. Can be gzipped.");
		parser.addArgument("-nr", "--nrandom")
			.type(Integer.class)
			.required(false)
			.setDefault(1)
			.help("Number of chains to generate. *Each* sequence in fasta"
					+ "input will be randomized this many times");
		parser.addArgument("-l", "--length")
			.type(Integer.class)
			.required(false)
			.setDefault(-1)
			.help("Lenght of each random sequence. Default to length of input sequence");

		parser.addArgument("-o", "--output")
			.type(String.class)
			.required(false)
			.setDefault("-")
			.help("Output file. Use - to stdout");
		parser.addArgument("-n", "--norder")
			.type(Integer.class)
			.required(false)
			.setDefault(1)
			.help("Order of the Markov model. Must be greater than 0.");
		parser.addArgument("-p", "--pct_pseudocounts")
			.type(Double.class)
			.required(false)
			.setDefault(0.1)
			.help("If the input sequence does not contain all the possible kmers,"
					+ "add the missing ones at this percentage of the total counts. "
					+ "E.g. -p 0.1 makes all the missing kmers to sum up to no more "
					+ "than 10% of the total actual counts. Setting to zero might terminate"
					+ "the chain with an error.");		
		parser.addArgument("--version", "-v").action(Arguments.version());
		Namespace opts= null;
		try{
			opts= parser.parseArgs(args);
		}
		catch(ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}		
		return(opts);
	}
}