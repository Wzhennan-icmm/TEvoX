#include "tevox.h"

#include <stdarg.h>
#include <stdio.h>

void tv_print_error(const char *format, ...)
{
    va_list args;

    fprintf(stderr, "tevox: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fputc('\n', stderr);
}

const char *tv_state_name(TvState state)
{
    static const char *names[] = {
        "PRESENT_ANNOTATED",
        "PRESENT_UNANNOTATED",
        "EMPTY_SITE_CONFIRMED",
        "STRUCTURAL_ALTERNATIVE",
        "FAMILY_OR_BOUNDARY_DISCORDANCE",
        "ASSEMBLY_GAP",
        "PROJECTION_AMBIGUOUS",
        "UNCALLABLE"
    };

    if (state < TV_PRESENT_ANNOTATED || state > TV_UNCALLABLE) {
        return "UNCALLABLE";
    }
    return names[state];
}

const char *tv_technical_state_name(TvTechnicalState state)
{
    static const char *names[] = {
        "CALLABLE", "GAP", "AMBIGUOUS", "UNCALLABLE"
    };

    if (state < TV_TECH_CALLABLE || state > TV_TECH_UNCALLABLE) {
        return "UNCALLABLE";
    }
    return names[state];
}

const char *tv_biological_state_name(TvBiologicalState state)
{
    static const char *names[] = {
        "PRESENT", "EMPTY", "STRUCTURAL_ALTERNATIVE", "UNKNOWN"
    };

    if (state < TV_BIO_PRESENT || state > TV_BIO_UNKNOWN) {
        return "UNKNOWN";
    }
    return names[state];
}

const char *tv_annotation_state_name(TvAnnotationState state)
{
    static const char *names[] = {
        "MATCHED", "MISSING", "FAMILY_CONFLICT", "NOT_APPLICABLE", "UNKNOWN"
    };

    if (state < TV_ANN_MATCHED || state > TV_ANN_UNKNOWN) {
        return "UNKNOWN";
    }
    return names[state];
}

const char *tv_claim_type_name(TvClaimType type)
{
    static const char *names[] = {
        "PRESENCE",
        "EMPTY_SITE",
        "STRUCTURAL_ALTERNATIVE",
        "ANNOTATION_DISCORDANCE",
        "NONE"
    };

    if (type < TV_CLAIM_PRESENCE || type > TV_CLAIM_NONE) {
        return "NONE";
    }
    return names[type];
}

const char *tv_identity_method_name(TvIdentityMethod method)
{
    static const char *names[] = {"MISSING", "EQX", "CS"};

    if (method < TV_IDENTITY_MISSING || method > TV_IDENTITY_CS) {
        return "MISSING";
    }
    return names[method];
}

const char *tv_evidence_origin_name(TvEvidenceOrigin origin)
{
    return origin == TV_EVIDENCE_NATIVE ? "NATIVE" : "DERIVED_REVERSE";
}

const char *tv_family_relation_name(TvFamilyRelation relation)
{
    static const char *names[] = {"MATCH", "CONFLICT", "UNKNOWN"};

    if (relation < TV_FAMILY_MATCH || relation > TV_FAMILY_UNKNOWN) {
        return "UNKNOWN";
    }
    return names[relation];
}
