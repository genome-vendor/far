 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_JOURNAL_GENERATED_FORWARDS_H 
#define SEQAN_HEADER_SEQUENCE_JOURNAL_GENERATED_FORWARDS_H 

//////////////////////////////////////////////////////////////////////////////
// NOTE: This file is automatically generated by build_forwards.py
//       Do not edit this file manually!
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// IndexOperatorValue

template <typename TValue, typename TSloppySpec> struct IndexOperatorValue;       	// "projects/library/seqan/sequence_journal/string_journal.h"(11)

//____________________________________________________________________________
// Journal

template <typename TValue, typename TSpec , typename TStringSpec , typename TSloppySpec > class Journal;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(39)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > class Journal;       	// "projects/library/seqan/sequence_journal/journal.h"(14)

//____________________________________________________________________________
// Sloppy

struct Sloppy;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(9)

//____________________________________________________________________________
// SloppyValue

template <typename TValue, typename TSloppySpec > struct SloppyValue;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(14)

//____________________________________________________________________________
// Strict

struct Strict;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(10)

//____________________________________________________________________________
// inorder_dbg_print

struct inorder_dbg_print;       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(130)

//____________________________________________________________________________
// inorder_offset

struct inorder_offset;       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(135)

//____________________________________________________________________________
// inorder_print

struct inorder_print;       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(125)

//____________________________________________________________________________
// jiter

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > class jiter;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(45)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > class jiter;       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(7)

//____________________________________________________________________________
// journal_iterator_proxy

template <typename TValue, typename TSloppySpec > struct journal_iterator_proxy;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(27)

//____________________________________________________________________________
// suffix_compare_functor

template <typename TString > struct suffix_compare_functor;       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(4)

//____________________________________________________________________________
// tree_visitor

struct tree_visitor;       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(114)

} //namespace seqan

//____________________________________________________________________________
// Deletion

struct Deletion;       	// "projects/library/seqan/sequence_journal/string_journal_operation.h"(65)

//____________________________________________________________________________
// Entry

struct Entry;       	// "projects/library/seqan/sequence_journal/string_journal_operation.h"(84)

//____________________________________________________________________________
// Insertion

struct Insertion;       	// "projects/library/seqan/sequence_journal/string_journal_operation.h"(46)

//____________________________________________________________________________
// JournalTest

template <typename TString > class JournalTest;       	// "projects/library/seqan/sequence_journal/string_journal_test_foundry.h"(28)

//____________________________________________________________________________
// Node

struct Node;       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(7)
struct Node;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(4)

//____________________________________________________________________________
// Operation

struct Operation;       	// "projects/library/seqan/sequence_journal/string_journal_operation.h"(4)

//____________________________________________________________________________
// tree_visitor

struct tree_visitor;       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(5)


//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// _lcp_length

template <typename TIteratorA, typename TIteratorB > inline size_t _lcp_length( TIteratorA it_a, TIteratorB it_b );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(30)

//____________________________________________________________________________
// _setLength

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void _setLength( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > * me, size_t new_length );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(287)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void _setLength( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, size_t new_length );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(294)

//____________________________________________________________________________
// _suffix_bigger

template <typename TIteratorA, typename TIteratorB > inline bool _suffix_bigger( TIteratorA it_a, TIteratorB it_b );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(18)

//____________________________________________________________________________
// append

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSource > void append( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, TSource const& source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(131)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSource > void append( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, TSource & source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(136)

//____________________________________________________________________________
// assign

template <typename TValue, typename TSpec, typename TStringSpec, typename TSourceSpec, typename TSloppySpec > void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, TSourceSpec > &source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(215)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(220)
template <typename TValue, typename TSpec, typename TStringSpec, typename TTargetSpec, typename TSloppySpec > void assign( String< TValue, TTargetSpec > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(225)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSourceSpec, typename TSloppySpec > void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &target, String< TValue, TSourceSpec > const &source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(230)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void assign( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const & source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(235)
template <typename TValue, typename TSpec, typename TStringSpec, typename TTargetSpec, typename TSloppySpec > void assign( String< TValue, TTargetSpec > &target, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &source );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(240)

//____________________________________________________________________________
// begin

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type begin( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, Standard );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(256)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type begin( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me, Standard );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(266)

//____________________________________________________________________________
// deleted

bool deleted( String< Pair< int, int > > const & deletions, int position );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(72)

//____________________________________________________________________________
// deleted_b

template <typename TString, typename TPos > inline bool deleted_b( TString const & deletions, TPos position );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(82)

//____________________________________________________________________________
// end

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type end( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, Standard );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(261)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > typename Iterator< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > >, Standard >::Type end( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me, Standard );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(271)

//____________________________________________________________________________
// erase

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition > void erase( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TPosition position );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(151)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition > void erase( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TPosition position, TPosition position_end );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(156)

//____________________________________________________________________________
// flatten

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void flatten( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(161)

//____________________________________________________________________________
// generate_shifts_and_deletions

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TIndex, typename TPos > inline void generate_shifts_and_deletions( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & string, String< TIndex > & indices, String< Pair< TPos, TPos > > & shifts, String< Pair< TPos, TPos > > & deletions );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(124)

//____________________________________________________________________________
// getValue

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPos > inline typename IndexOperatorValue< TValue, TSloppySpec >::Type getValue( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, TPos pos );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(200)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPos > inline TValue const & getValue( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const & me, TPos pos );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(205)

//____________________________________________________________________________
// get_shift

inline int get_shift( String< Pair< int, int > > const & limits, int position );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(43)

//____________________________________________________________________________
// get_shift_b

template <typename TString, typename TPos > inline int get_shift_b( TString const & limits, TPos position );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(56)

//____________________________________________________________________________
// id

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void const * id( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(281)

//____________________________________________________________________________
// insert

template <typename TValue, typename TSpec, typename TStringSpec, typename TString, typename TSloppySpec > void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TString &insert_string );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(113)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TValue value );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(118)
template <typename TValue, typename TSpec, typename TStringSpec, typename TIterator, typename TSloppySpec > void insert( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TIterator it, size_t number );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(126)

//____________________________________________________________________________
// j_goDownLeft

template <typename TIter > inline bool j_goDownLeft( TIter &it );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(306)

//____________________________________________________________________________
// j_goDownRight

template <typename TIter > inline bool j_goDownRight( TIter &it );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(291)

//____________________________________________________________________________
// j_goNext

template <typename TIter > inline bool j_goNext( TIter &it );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(273)

//____________________________________________________________________________
// j_goPrev

template <typename TIter > inline bool j_goPrev( TIter &it );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(343)

//____________________________________________________________________________
// j_goUp

template <typename TIter > inline size_t j_goUp( TIter &it );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(321)

//____________________________________________________________________________
// length

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > size_t length( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > const &me );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(276)

//____________________________________________________________________________
// operator*

template <typename TValue, typename TSloppySpec > typename SloppyValue< TValue, TSloppySpec >::Type operator*( journal_iterator_proxy< TValue, TSloppySpec > jtp );       	// "projects/library/seqan/sequence_journal/string_journal_forwards.h"(34)

//____________________________________________________________________________
// printSA

template <typename TIndex, typename TString > void printSA( TIndex & index, TString & string );       	// "projects/library/seqan/sequence_journal/string_journal_debug.h"(3)

//____________________________________________________________________________
// printSA_LCP

template <typename TIndex, typename TString > void printSA_LCP( TIndex & index, TString & string );       	// "projects/library/seqan/sequence_journal/string_journal_debug.h"(14)

//____________________________________________________________________________
// print_pairs

template <typename TValue1, typename TValue2, typename TTag > void print_pairs( String< Pair< TValue1, TValue2, TTag > > const & string );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(116)

//____________________________________________________________________________
// replace

template <typename TValue, typename TSpec, typename TStringSpec, typename TString, typename TSloppySpec > void replace( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TString &insert_string );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(141)
template <typename TValue, typename TSpec, typename TStringSpec, typename TIterator, typename TSloppySpec > void replace( size_t position, String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &journal_string, TIterator it, size_t number );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(146)

//____________________________________________________________________________
// reserve

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TSize, typename TExpand> inline typename Size< String<TValue, Journal < TSpec, TStringSpec, TSloppySpec > > >::Type reserve( String<TValue, Journal < TSpec, TStringSpec, TSloppySpec > > & me, TSize new_capacity, Tag<TExpand> const tag);       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(210)

//____________________________________________________________________________
// resize

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TLength, typename TExpand > inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type resize( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, TLength new_length, Tag<TExpand> const );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(189)

//____________________________________________________________________________
// resizeSpace

template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec> inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type resizeSpace( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > &me, size_t size, size_t pos_begin, size_t pos_end, size_t limit );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(171)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec, typename TPosition, typename TExpand> inline typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type resizeSpace( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & me, typename Size< String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > >::Type size, TPosition pos_begin, TPosition pos_end, Tag<TExpand> const);       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(183)

//____________________________________________________________________________
// synchronize_ChildTab

template <typename TIndex, typename TString > void synchronize_ChildTab( TIndex & index, TString & string, std::vector<Entry>& entries );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(753)

//____________________________________________________________________________
// synchronize_SA_LCP

template <typename TIndex, typename TString > void synchronize_SA_LCP( TIndex & index, TString & string, std::vector<Entry>& entries );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(621)

//____________________________________________________________________________
// synchronize_index

template <typename TIndex, typename TString, typename TSpec > inline void synchronize_index( TIndex & index, StringSet< TString, TSpec > & stringset );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(159)
template <typename TIndex, typename TString > inline void synchronize_index( TIndex & index, TString & string );       	// "projects/library/seqan/sequence_journal/string_journal_utility.h"(466)

//____________________________________________________________________________
// tree_visitor::~tree_visitor

tree_visitor::~tree_visitor();       	// "projects/library/seqan/sequence_journal/string_journal_node.h"(123)

//____________________________________________________________________________
// value

template <typename TValue, typename TSpec, typename TStringSpec, typename TPosition > TValue const & value( String< TValue, Journal< TValue, TSpec, TStringSpec, Strict > > const & string, TPosition position );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(19)
template <typename TValue, typename TSpec, typename TStringSpec, typename TPosition > TValue & value( String< TValue, Journal< TValue, TSpec, TStringSpec, Sloppy > > & string, TPosition position );       	// "projects/library/seqan/sequence_journal/string_journal_interface.h"(24)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > inline TValue const & value( jiter< TValue, TSpec, TStringSpec, TSloppySpec > const & me );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(263)
template <typename TValue, typename TSpec, typename TStringSpec, typename TSloppySpec > inline TValue const & value( jiter< TValue, TSpec, TStringSpec, TSloppySpec > & me );       	// "projects/library/seqan/sequence_journal/iterator_journal.h"(268)

} //namespace seqan

//____________________________________________________________________________
// generate_random_string

template <typename TSpec > void generate_random_string( unsigned int length, seqan::String< char, TSpec > & the_string, seqan::String< char > const & alphabet_string );       	// "projects/library/seqan/sequence_journal/string_journal_test_foundry.h"(14)

//____________________________________________________________________________
// randomNumber

unsigned int randomNumber( int upper_bound );       	// "projects/library/seqan/sequence_journal/string_journal_test_foundry.h"(4)

#endif

