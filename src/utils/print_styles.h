#pragma once

#include <string>

#include "glog/logging.h"

namespace circinus {

// font style
#define BOLD(str) "\e[1m" << str << "\e[0m"
#define UNDERLINE(str) "\e[4m" << str << "\e[0m"
#define BLINK(str) "\e[5m" << str << "\e[0m"
#define INVERSE(str) "\e[7m" << str << "\e[0m"
// font color
#define RED(str) "\e[31m" << str << "\e[0m"
#define GREEN(str) "\e[32m" << str << "\e[0m"
#define YELLOW(str) "\e[33m" << str << "\e[0m"
#define BLUE(str) "\e[34m" << str << "\e[0m"
#define MAGENTA(str) "\e[35m" << str << "\e[0m"
#define CYAN(str) "\e[36m" << str << "\e[0m"
// background color
#define BRED(str) "\e[41m" << str << "\e[0m"
#define BGREEN(str) "\e[42m" << str << "\e[0m"
#define BYELLOW(str) "\e[43m" << str << "\e[0m"
#define BBLUE(str) "\e[44m" << str << "\e[0m"
#define BMAGENTA(str) "\e[45m" << str << "\e[0m"
#define BCYAN(str) "\e[46m" << str << "\e[0m"

// font style
#define LogBold(code) LOG(INFO) << "\e[1m" << code << "\e[0m"
#define LogUnderline(code) LOG(INFO) << "\e[4m" << code << "\e[0m"
#define LogBlink(code) LOG(INFO) << "\e[5m" << code << "\e[0m"
#define LogInverse(code) LOG(INFO) << "\e[7m" << code << "\e[0m"
// font color
#define LogRed(code) LOG(INFO) << "\e[31m" << code << "\e[0m"
#define LogGreen(code) LOG(INFO) << "\e[32m" << code << "\e[0m"
#define LogYellow(code) LOG(INFO) << "\e[33m" << code << "\e[0m"
#define LogBlue(code) LOG(INFO) << "\e[34m" << code << "\e[0m"
#define LogMagenta(code) LOG(INFO) << "\e[35m" << code << "\e[0m"
#define LogCyan(code) LOG(INFO) << "\e[36m" << code << "\e[0m"
// background color
#define LogRedBackground(code) LOG(INFO) << "\e[41m" << code << "\e[0m"
#define LogGreenBackground(code) LOG(INFO) << "\e[42m" << code << "\e[0m"
#define LogYellowBackground(code) LOG(INFO) << "\e[43m" << code << "\e[0m"
#define LogBlueBackground(code) LOG(INFO) << "\e[44m" << code << "\e[0m"
#define LogMagentaBackground(code) LOG(INFO) << "\e[45m" << code << "\e[0m"
#define LogCyanBackground(code) LOG(INFO) << "\e[46m" << code << "\e[0m"

}  // namespace circinus
