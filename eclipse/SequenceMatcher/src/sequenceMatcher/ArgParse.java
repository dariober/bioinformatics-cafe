package sequenceMatcher;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;


public class ArgParse {
	
	public static String VERSION= "0.1.0";
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		
		String mainHelp= "DESCRIPTION\n"
				+ "This program answers the question: Which sequences in file a "
				+ "are similar to the sequences in file b.\n"
				+ "\n"
				+ "Output format:\n"
				+ "1. Sequence name from the A file\n"
				+ "2. Sequence name from the B file\n"
				+ "3. Strand (+: Sequence matched as they are; -: Sequence A reverse complemented)\n"
				+ "4. Edit distance\n"
				+ "5. Length sequence A\n"
				+ "6. Length sequence B\n"
				+ "7. Sequence A, if matched as rev. comp. it is rev comp'd here.\n"
				+ "8. Sequence B\n"
				+ "\n"
				+ "NB: Names of fasta sequences must be unique within file.";
		
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("SequenceMatcher")
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description(mainHelp);

		parser.addArgument("--version", "-v").action(Arguments.version());

		Subparsers subparsers= parser.addSubparsers()
				.dest("subcmd")
				.description("Matching & parsing")
				.help("Run SequenceMatcher match -h to view options.");
		
		/* --------------------------------------------------------------------- */

		Subparser matchSubparser= subparsers.addParser("match");		
		matchSubparser.description("Main program to perform matching and alignment.");
		
		matchSubparser.addArgument("-a", "--a")
			.type(String.class)
			.required(true)
			.help("Fasta file 'A'; can be gzip'd. Use '-' to read from stdin. "
					+ "In SAM format this is 'Reference'");
		
		matchSubparser.addArgument("-b", "--b")
			.type(String.class)
			.required(true)
			.help("Fasta file 'B'; can be gzip'd. Use '-' to read from stdin. "
					+ "In SAM format these are 'Reads'");
		
		matchSubparser.addArgument("-m", "--method")
			.type(String.class)
			.required(false)
			.setDefault("LD")
			.choices("LD", "HD")
			.help("Method to determine the threshold edit distance: LD (Levenshtein) or HD (Hamming)"
					+ "Hamming is much faster.");
		
		matchSubparser.addArgument("-nm", "--nm")
			.type(Integer.class)
			.setDefault(-1)
			.help("Maximum edit distance to output a match. If -1 (default) all sequence pairs are returned.");
		
		matchSubparser.addArgument("-norc", "--norc")
			.action(Arguments.storeTrue())
			.help("Do not reverse complement the sequences in file A. I.e. only match sequences as they are.");	

		matchSubparser.addArgument("-noaln", "--noaln")
			.action(Arguments.storeTrue())
			.help("Do not align sequence, just match them (faster).");	

		matchSubparser.addArgument("-noLD", "--noLD")
			.action(Arguments.storeTrue())
			.help("Do not compute Levenshtein distance (faster).");	

		matchSubparser.addArgument("-noJWD", "--noJWD")
			.action(Arguments.storeTrue())
			.help("Do not compute Jaro-Winkler distance (faster).");	

		matchSubparser.addArgument("-of", "--outfmt")
			.type(String.class)
			.required(false)
			.setDefault("tab")
			.choices("tab", "sam")
			.help("Output format. 'tab': Tab delim (see above); 'sam' SAM format");	

		/* -------------------------------------------------------------------- */
		
		Subparser convertSubparser= subparsers.addParser("convert");
		convertSubparser.description("Convert tab to sam and viceversa");
		
		convertSubparser.addArgument("-i", "--input")
			.type(String.class)
			.required(true)
			.help("Input file to convert. can be gzip'd. Use '-' to read from stdin.");
		
		convertSubparser.addArgument("-of", "--outfmt")
			.type(String.class)
			.choices("tab", "sam")
			.required(true)
			.help("Output format. 'tab' to convert sam input to tab; 'sam' viceversa.");

		convertSubparser.addArgument("-a", "--a")
		.type(String.class)
		.help("Required for converting tab to sam: Reference fasta file aka as 'a' file; "
				+ "can be gzip'd. Use '-' to read from stdin. ");

		/* --------------------------------------------------------------------- */
		
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
	
	public static void validateArgs(Namespace opts){
		//		
	}
}
