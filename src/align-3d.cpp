// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// 3d dingen

#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>
#include <numeric>
#include <vector>
#include <map>
#include <valarray>

#include <cif++/Cif++.hpp>

// #include <sys/times.h>
// #include <sys/resource.h>

#include "align-3d.hpp"
#include "matrix.hpp"

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

// namespace fs = boost::filesystem;
// namespace io = boost::iostreams;
namespace ba = boost::algorithm;
// namespace bm = boost::math;

using mmcif::Point;
using mmcif::quaternion;
using mmcif::Centroid;
using mmcif::QuaternionToAngleAxis;

// --------------------------------------------------------------------

double CalculateRMSD(const std::vector<mmcif::Point>& pa, const std::vector<mmcif::Point>& pb)
{
	return RMSd(pa, pb);
}

void AlignIterative(std::vector<Point>& cAlphaA, std::vector<Point>& cAlphaB,
	Point& outTranslationA, Point& outTranslationB, quaternion& outRotation)
{
	if (cAlphaA.size() != cAlphaB.size())
		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
	if (cAlphaA.size() < 3)
		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");

	outTranslationA = Centroid(cAlphaA);
	for (Point& pt: cAlphaA)
		pt -= outTranslationA;
	
	outTranslationB = Centroid(cAlphaB);
	for (Point& pt: cAlphaB)
		pt -= outTranslationB;
	
	outRotation = AlignPoints(cAlphaA, cAlphaB);
	for (Point& pt: cAlphaB)
		pt.rotate(outRotation);

	double angle; Point axis;
	tie(angle, axis) = QuaternionToAngleAxis(outRotation);
	if (cif::VERBOSE > 1)
		cerr << "before iterating, rotation: " << angle << " degrees rotation around axis " << axis << endl;
}

// bool align_proteins2(MProtein& a, MProtein& b, char chainA, char chainB,
// 	Point& outTranslationA, Point& outTranslationB, MQuaternion& outRotation)
// {
// 	sequence sa, sb;
	
// 	a.GetSequence(chainA, sa);
// 	b.GetSequence(chainB, sb);

// 	vector<Point> cAlphaA, cAlphaB;
	
// 	a.GetCAlphaLocations(chainA, cAlphaA);
// 	b.GetCAlphaLocations(chainB, cAlphaB);

// 	uint32 N = cAlphaA.size();
// 	uint32 M = cAlphaB.size();
	
// 	matrix<int32> B(N, M, 0);
// 	matrix<int8> traceback(N, M, 2);
// 	int32 high = 0, highA, highB;
	
// 	for (uint32 ai = 0; ai < cAlphaA.size(); ++ai)
// 	{
// 		for (uint32 bi = 0; bi < cAlphaB.size(); ++bi)
// 		{
// 			float d = Distance(cAlphaA[ai], cAlphaB[bi]);

// 			int32 v = 0;
// 			if (d < 3.5)
// 				v = 1;

// 			if (ai > 0 and bi > 0)
// 				v += B(ai - 1, bi - 1);
			
// 			int32 ga = -1;
// 			if (ai > 0)
// 				ga = B(ai - 1, bi) - 1;
			
// 			int32 gb = -1;
// 			if (bi > 0)
// 				gb = B(ai, bi - 1) - 1;
			
// 			if (VERBOSE > 5)
// 				cerr << "d(" << ai << ',' << bi << ") = " << d << " ga: " << ga << " gb: " << gb << " v: " << v << endl;

// 			if (v >= ga and v >= gb)
// 			{
// 				B(ai, bi) = v;
// 				traceback(ai, bi) = 0;
// 			}
// 			else if (ga > gb)
// 			{
// 				B(ai, bi) = ga;
// 				traceback(ai, bi) = 1;
// 			}
// 			else
// 			{
// 				B(ai, bi) = gb;
// 				traceback(ai, bi) = -1;
// 			}
			
// 			if (B(ai, bi) > high)
// 			{
// 				high = B(ai, bi);
// 				highA = ai;
// 				highB = bi;
// 			}
// 		}
// 	}

// 	if (VERBOSE > 4)
// 		print_matrix(cerr, traceback, sa, sb);
	
// 	vector<Point> newCAlphaA, newCAlphaB;
	
// 	int32 ai = highA, bi = highB;
	
// 	while (ai > 0 and bi > 0)
// 	{
// 		switch (traceback(ai, bi))
// 		{
// 			case 0:
// 			{
// 				float d = Distance(cAlphaA[ai], cAlphaB[bi]);
// 				if (d < 3.5)
// 				{
// 					if (VERBOSE > 3)
// 						cerr << "map " << ai << '(' << kAA[sa[ai]]
// 							 << ") to " << bi << '(' << kAA[sb[bi]] << ')' << endl;
// 					newCAlphaA.push_back(cAlphaA[ai]);
// 					newCAlphaB.push_back(cAlphaB[bi]);
// 				}
// 				--ai;
// 				--bi;
// 				break;
// 			}
			
// 			case 1:
// 				--ai;
// 				break;
			
// 			case -1:
// 				--bi;
// 				break;
			
// 			default:
// 				assert(false);
// 				break;
// 		}
// 	}

// 	outTranslationA = Centroid(newCAlphaA);
// 	foreach (Point& pt, newCAlphaA)
// 		pt -= outTranslationA;
	
// 	outTranslationB = Centroid(newCAlphaB);
// 	foreach (Point& pt, newCAlphaB)
// 		pt -= outTranslationB;
	
// 	outRotation = AlignPoints(newCAlphaA, newCAlphaB);

// 	double angle;
// 	Point axis;

// 	tr1::tie(angle, axis) = QuaternionToAngleAxis(outRotation);
// 	if (VERBOSE > 1)
// 		cerr << "  translation: " << outTranslationA << " and " << outTranslationB << endl
// 			 << "  rotation: " << angle << " degrees around axis " << axis << endl;

// 	return angle > 0.01;
// }

// void align_proteins(MProtein& a, char chainA, MProtein& b, char chainB,
// 	uint32 iterations,
// 	substitution_matrix_family& mat, float gop, float gep, float magic)
// {
// 	vector<Point> cAlphaA, cAlphaB;

// 	a.CalculateSecondaryStructure();
// 	b.CalculateSecondaryStructure();

// 	// fetch sequences
// 	entry ea(1, a.GetID()), eb(2, b.GetID());
// 	a.GetSequence(chainA, ea);
// 	b.GetSequence(chainB, eb);
	
// 	vector<entry*> aa; aa.push_back(&ea);
// 	vector<entry*> ab; ab.push_back(&eb);

// 	vector<entry*> ac;
// 	joined_node n(new leaf_node(ea), new leaf_node(eb), 0.1, 0.1);

// 	align(&n, aa, ab, ac, mat, gop, gep, magic, true);

// 	// now based on this alignment, select c-alpha's that we try to align
// 	assert(ac.size() == 2);
// 	assert(ac.front()->m_seq.length() == ac.back()->m_seq.length());
	
// 	const sequence& sa = ac.front()->m_seq;
// 	const sequence& sb = ac.back()->m_seq;
	
// 	const vector<int16>& pa = ac.front()->m_positions;
// 	const vector<int16>& pb = ac.back()->m_positions;

// 	sequence nsa, nsb;
	
// 	for (uint32 i = 0; i < sa.length(); ++i)
// 	{
// 		if (sa[i] == kSignalGapCode or sb[i] == kSignalGapCode)
// 			continue;

// 		nsa.push_back(sa[i]);
// 		cAlphaA.push_back(a.GetCAlphaPosition(chainA, pa[i]));

// 		nsb.push_back(sb[i]);
// 		cAlphaB.push_back(b.GetCAlphaPosition(chainB, pb[i]));
// 	}
	
// 	if (cAlphaA.size() != cAlphaB.size())
// 		throw logic_error("Protein A and B should have the same number of c-alpha atoms");
	
// 	if (cAlphaA.size() < 3)
// 		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");
	
// 	Point ta, tb;
// 	MQuaternion rotation;
	
// 	// 
	
// 	if (cAlphaA.size() != cAlphaB.size())
// 		throw logic_error("Protein A and B should have the same number of c-alhpa atoms");
	
// 	if (cAlphaA.size() < 3)
// 		throw logic_error("Protein A and B should have at least 3 of c-alpha atoms");

// 	//

// 	ta = Centroid(cAlphaA);
// 	foreach (Point& pt, cAlphaA)
// 		pt -= ta;
	
// 	tb = Centroid(cAlphaB);
// 	foreach (Point& pt, cAlphaB)
// 		pt -= tb;
	
// 	rotation = AlignPoints(cAlphaA, cAlphaB);
// 	foreach (Point& pt, cAlphaB)
// 		pt.Rotate(rotation);

// 	double angle; Point axis;
// 	tr1::tie(angle, axis) = QuaternionToAngleAxis(rotation);
// 	if (VERBOSE > 1)
// 		cerr << "before iterating, rotation: " << angle << " degrees rotation around axis " << axis << endl;
	
// 	uint32 iteration = 0;
// 	for (;;)
// 	{
// 		a.Translate(-ta);
// 		b.Translate(-tb);
// 		b.Rotate(rotation);
		
// 		if (VERBOSE > 1)
// 			cerr << "RMSd: (" << iteration << ") " << CalculateRMSD(a, b, chainA, chainB) << endl;

// 		if (iteration++ >= iterations)
// 			break;
		
// 		if (not align_proteins2(a, b, chainA, chainB, ta, tb, rotation))
// 			break;
// 	}
// }

// void create_entries(MProtein& a, MProtein& b, char chainA, char chainB,
// 	entry& ea, entry& eb)
// {
// 	a.GetSequence(chainA, ea);
// 	b.GetSequence(chainB, eb);
	
// 	ea.m_positions = vector<int16>(ea.m_seq.length());
// 	eb.m_positions = vector<int16>(eb.m_seq.length());

// 	vector<Point> cAlphaA, cAlphaB;
	
// 	a.GetCAlphaLocations(chainA, cAlphaA);
// 	b.GetCAlphaLocations(chainB, cAlphaB);

// 	uint32 N = cAlphaA.size();
// 	uint32 M = cAlphaB.size();
	
// 	matrix<int32> B(N, M, 0);
// 	matrix<int8> traceback(N, M, 2);
// 	int32 high = 0, highA, highB;
	
// 	for (uint32 ai = 0; ai < cAlphaA.size(); ++ai)
// 	{
// 		for (uint32 bi = 0; bi < cAlphaB.size(); ++bi)
// 		{
// 			float d = Distance(cAlphaA[ai], cAlphaB[bi]);

// 			int32 v = 0;
// 			if (d < 3.5)
// 				v = 1;

// 			if (ai > 0 and bi > 0)
// 				v += B(ai - 1, bi - 1);
			
// 			int32 ga = -1;
// 			if (ai > 0)
// 				ga = B(ai - 1, bi) - 1;
			
// 			int32 gb = -1;
// 			if (bi > 0)
// 				gb = B(ai, bi - 1) - 1;
			
// 			if (v >= ga and v >= gb)
// 			{
// 				B(ai, bi) = v;
// 				traceback(ai, bi) = 0;
// 			}
// 			else if (ga > gb)
// 			{
// 				B(ai, bi) = ga;
// 				traceback(ai, bi) = 1;
// 			}
// 			else
// 			{
// 				B(ai, bi) = gb;
// 				traceback(ai, bi) = -1;
// 			}
			
// 			if (B(ai, bi) > high)
// 			{
// 				high = B(ai, bi);
// 				highA = ai;
// 				highB = bi;
// 			}
// 		}
// 	}

// 	int32 ai = highA, bi = highB;
// 	int16 posNr = max(N, M) + 1;
	
// 	while (ai > 0 and bi > 0)
// 	{
// 		switch (traceback(ai, bi))
// 		{
// 			case 0:
// 				if (Distance(cAlphaA[ai], cAlphaB[bi]) <= 1.9)
// 				{
// 					assert(posNr > 0);
// 					ea.m_positions[ai] = posNr;
// 					eb.m_positions[bi] = posNr;
// 					--posNr;
// 				}
// 				--ai;
// 				--bi;
// 				break;
			
// 			case 1:
// 				--ai;
// 				break;
			
// 			case -1:
// 				--bi;
// 				break;
			
// 			default:
// 				assert(false);
// 				break;
// 		}
// 	}
// }

// void align_structures(
// 	std::istream& structureA, std::istream& structureB,
// 	char chainA, char chainB,
// 	uint32 iterations,
// 	substitution_matrix_family& mat, float gop, float gep, float magic,
// 	std::vector<entry*>& alignment)
// {
// //	string pdbid_a = structureA;
// //	if (pdbid_a.length() == 5)
// //	{
// //		chainA = pdbid_a[4];
// //		pdbid_a.erase(4, 5);
// //	}
// //	else if (pdbid_a.length() != 4)
// //		throw mas_exception("Please specify a PDB ID in 4 letter code");
// //
// //	string pdbid_b = structureB;
// //	if (pdbid_b.length() == 5)
// //	{
// //		chainB = pdbid_b[4];
// //		pdbid_b.erase(4, 5);
// //	}
// //	else if (pdbid_b.length() != 4)
// //		throw mas_exception("Please specify a PDB ID in 4 letter code");

// //	stringstream file_a(pdb->GetDocument(pdbid_a));
// 	MProtein a(structureA, false);
	
// 	if (chainA == 0)
// 		chainA = a.GetFirstChainID();

// //	stringstream file_b(pdb->GetDocument(pdbid_b));
// 	MProtein b(structureB, false);

// 	if (chainB == 0)
// 		chainB = b.GetFirstChainID();

// 	align_proteins(a, chainA, b, chainB, iterations, mat, gop, gep, magic);
	
// 	if (VERBOSE)
// 		cerr << "RMSD: " << CalculateRMSD(a, b, chainA, chainB) << endl;

// 	entry* ea = new entry(1, a.GetID());
// 	entry* eb = new entry(2, b.GetID());
// 	create_entries(a, b, chainA, chainB, *ea, *eb);
	
// 	vector<entry*> aa, ab;
// 	aa.push_back(ea);
// 	ab.push_back(eb);
// 	joined_node n(new leaf_node(*ea), new leaf_node(*eb), 0.1, 0.1);

// 	align(&n, aa, ab, alignment, mat, gop, gep, magic, false);
	
// 	MProtein c;

// 	c.SetChain('A', a.GetChain(chainA));
// 	c.SetChain('B', b.GetChain(chainB));
	
// 	string pdbid_a = a.GetID();
// 	if (isprint(chainA) and not isspace(chainA))
// 		pdbid_a += chainA;
	
// 	string pdbid_b = b.GetID();
// 	if (isprint(chainB) and not isspace(chainB))
// 		pdbid_b += chainB;
	
// 	ofstream file_o((pdbid_a + '-' + pdbid_b + ".pdb").c_str());
// 	c.WritePDB(file_o);
// }

