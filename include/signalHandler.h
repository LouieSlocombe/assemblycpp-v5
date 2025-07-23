/**
 * @file signalHandler.h
 * @brief ccode to handle interrupt signal and output best pathway found so far
 */
#pragma once
#ifdef _WIN32
#include <windows.h>
/**
 * @brief handles interrupt signal so that best pathway found can be output
 */
BOOL CtrlHandler(DWORD fdwCtrlType);
#else
#include <csignal>
/**
 * @brief handles interrupt signal so that best pathway found can be output
 */
void signalHandler(int signum);
#endif