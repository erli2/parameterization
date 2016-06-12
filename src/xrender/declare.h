#pragma once

#ifdef XRENDER_PROJ
#define DEVICE_API __declspec(dllexport)
#else
#define DEVICE_API __declspec(dllimport)
#endif


