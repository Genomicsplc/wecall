// All content Copyright (C) 2018 Genomics plc
#ifndef TIMER_HPP
#define TIMER_HPP
#include <chrono>
#include <string>
#include <map>
#include <memory>

namespace echidna
{
namespace utils
{
    std::string encodeString( std::string rawString );

    class Timer final
    {
    public:
        Timer( std::string type, std::map< std::string, std::string > metadata );
        ~Timer();

        void start();
        void pause();

    private:
        std::string format_metadata() const;
        std::string m_type;
        std::map< std::string, std::string > m_metadata;

        std::chrono::steady_clock::time_point m_start;
        long m_duration;
    };

    using timerPtr_t = std::shared_ptr< Timer >;

    std::map< std::string, std::string > fileMetaData( std::string filename );

    class ScopedTimerTrigger
    {
    public:
        ScopedTimerTrigger( timerPtr_t timer ) : m_timer( timer ) { m_timer->start(); }

        ~ScopedTimerTrigger() { m_timer->pause(); }

    private:
        timerPtr_t m_timer;
    };
}
}

#endif
