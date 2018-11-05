// All content Copyright (C) 2018 Genomics plc
#ifndef ALIGN_SCORER_HPP
#define ALIGN_SCORER_HPP

namespace wecall
{
/// All alignment algorithms and associated classes live in
/// this namespace
namespace alignment
{
    ///----------------------------------------------------------------------------------------
    /// Functor-style class to score alignment matches and gaps.
    ///----------------------------------------------------------------------------------------
    class AlignScorer
    {
    public:
        AlignScorer() = delete;
        explicit AlignScorer( const int gapPenalty ) : m_gapPenalty( gapPenalty ) {}
        explicit AlignScorer( const AlignScorer & rhs ) = delete;

        int scoreMatch( const char x, const char y ) const { return x == y ? 1 : -1; }
        int scoreGap() const { return m_gapPenalty; }

    private:
        const int m_gapPenalty;  /// Gap extension penalty
    };
}
}

#endif
