//
// Created by Admin on 2019/1/11.
//

#ifndef HTS_OBOE_SINEGENERATOR_H
#define HTS_OBOE_SINEGENERATOR_H

#include <math.h>
#include <cstdint>

class SineGenerator
{
public:
    SineGenerator();
    ~SineGenerator() = default;

    void setup(double frequency, int32_t frameRate);

    void setup(double frequency, int32_t frameRate, float amplitude);

    void setSweep(double frequencyLow, double frequencyHigh, double seconds);

    void render(int16_t *buffer, int32_t channelStride, int32_t numFrames);

    void render(float *buffer, int32_t channelStride, int32_t numFrames);

private:
    double mAmplitude;
    double mPhase = 0.0;
    int32_t mFrameRate;
    double mPhaseIncrement;
    double mPhaseIncrementLow;
    double mPhaseIncrementHigh;
    double mUpScaler = 1.0;
    double mDownScaler = 1.0;
    bool   mGoingUp = false;
    bool   mSweeping = false;

    void advancePhase();

    double getPhaseIncremement(double frequency);
};

#endif //HTS_OBOE_SINEGENERATOR_H
