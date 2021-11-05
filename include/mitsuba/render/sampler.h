#pragma once

#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/random.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/mitsuba.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER Sampler : public Object {
public:
    MTS_IMPORT_TYPES()

    /**
     * \brief Create a clone of this sampler
     *
     * The clone is allowed to be different to some extent, e.g. a pseudorandom
     * generator should be based on a different random seed compared to the
     * original. All other parameters are copied exactly.
     *
     * May throw an exception if not supported. Cloning may also change the
     * state of the original sampler (e.g. by using the next 1D sample as a
     * seed for the clone).
     */
    virtual ref<Sampler> clone() = 0;

    /**
     * \brief Deterministically seed the underlying RNG, if applicable.
     *
     * In the context of wavefront ray tracing & dynamic arrays, this function
     * must be called with a \c seed_value matching the size of the wavefront.
     */
    virtual void seed(uint64_t seed_offset);

    /**
     * \brief Advance to the next sample.
     *
     * A subsequent call to \c next_1d or \c next_2d will access the first
     * 1D or 2D components of this sample.
     */
    virtual void advance();

    /// Retrieve the next component value from the current sample
    virtual Float next_1d(Mask active = true);

    /// Retrieve the next two component values from the current sample
    virtual Point2f next_2d(Mask active = true);

    /// Return the number of samples per pixel
    size_t sample_count() const { return m_sample_count; }

    /// Return the size of the wavefront (or 0, if not seeded)
    virtual size_t wavefront_size() const = 0;

    /// Return whether the sampler was seeded
    // bool seeded() const { return m_wavefront_size > 0; }

    /// Set the number of samples per pass in wavefront modes (default is 1)
    void set_samples_per_wavefront(uint32_t samples_per_wavefront);

    MTS_DECLARE_CLASS()
protected:
    Sampler(const Properties &props);
    virtual ~Sampler();

    /// Generates a array of seeds where the seed values are unique per sequence
    UInt32 compute_per_sequence_seed(uint32_t seed_offset) const;
    /// Return the index of the current sample
    UInt32 current_sample_index() const;

protected:
    /// Base seed value
    uint64_t m_base_seed;
    /// Number of samples per pixel
    uint32_t m_sample_count;
    /// Number of samples per pass in wavefront modes (default is 1)
    uint32_t m_samples_per_wavefront;
    /// Size of the wavefront (or 0, if not seeded)
    uint32_t m_wavefront_size;
    /// Index of the current dimension in the sample
    uint32_t m_dimension_index;
    /// Index of the current sample in the sequence
    uint32_t m_sample_index;

    // size_t m_sample_count;
    // ScalarUInt64 m_base_seed;
};

// /// Interface for sampler plugins based on the PCG32 random number generator
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER PCG32Sampler : public Sampler<Float, Spectrum> {
// public:
//     MTS_IMPORT_BASE(Sampler, m_base_seed)
//     MTS_IMPORT_TYPES()
//     using PCG32 = mitsuba::PCG32<UInt32>;

//     virtual void seed(uint64_t seed_offset) override;
//     virtual size_t wavefront_size() const override;

//     MTS_DECLARE_CLASS()
// protected:
//     PCG32Sampler(const Properties &props);

// protected:
//     std::unique_ptr<PCG32> m_rng;
// };

template <typename Float, typename Spectrum>
class LowDiscrepancySampler final : public Sampler<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Sampler, m_sample_count, m_base_seed, // seeded,
                    m_samples_per_wavefront, m_dimension_index,
                    current_sample_index, compute_per_sequence_seed)
    MTS_IMPORT_TYPES()

    LowDiscrepancySampler(const Properties &props = Properties());

    LowDiscrepancySampler(int res_s, const Properties &props = Properties());

    ref<Sampler<Float, Spectrum>> clone() override {
        LowDiscrepancySampler *sampler   = new LowDiscrepancySampler();
        sampler->m_sample_count          = m_sample_count;
        sampler->m_samples_per_wavefront = m_samples_per_wavefront;
        sampler->m_base_seed             = m_base_seed;
        return sampler;
    }

    // void seed(uint64_t seed_offset, size_t wavefront_size) override {
    //     Base::seed(seed_offset, wavefront_size);
    //     m_scramble_seed = compute_per_sequence_seed(seed_offset);
    // }

    void seed(uint64_t seed_offset) override;

    // TODO: I am not really sure what this does
    virtual size_t wavefront_size() const override { return 1; }

    Float next_1d_test();

    Point2f next_2d_test();

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "LowDiscrepancySampler [" << std::endl
            << "  sample_count = " << m_sample_count << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    /// Per-sequence scramble seed
    uint32_t m_scramble_seed;
};

MTS_EXTERN_CLASS_RENDER(Sampler)
MTS_EXTERN_CLASS_RENDER(LowDiscrepancySampler)
// MTS_EXTERN_CLASS_RENDER(PCG32Sampler)
NAMESPACE_END(mitsuba)
